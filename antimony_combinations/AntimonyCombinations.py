import site, os, glob
import pandas, numpy
import re
import tellurium as te

# site.addsitedir(r'/home/ncw135/Documents/pycotools3')
# site.addsitedir(r'D:\pycotools3')
from pycotools3 import model, tasks, viz
from itertools import combinations
from collections import OrderedDict
import matplotlib.pyplot as plt
import seaborn
import yaml
import logging
from copy import deepcopy

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


class HypothesisExtension:
    """
    Data class for storing information about a hypothesis extension. For usage
    see :py:class:`Combinations`.
    """
    def __init__(self, name, reaction,
                 rate_law, mode='additive',
                 to_replace=None):
        self.name = name
        self.reaction = reaction
        self.rate_law = rate_law
        self.mode = mode
        self.to_replace = to_replace

        for i in [self.name, self.reaction, self.rate_law, self.mode]:
            if not isinstance(i, str):
                raise ValueError('attribute "{}" should be a string, not {}'.format(i, type(i)))

    def __str__(self):
        return f'{self.name}: {self.reaction}; {self.rate_law}'

    def __repr__(self):
        return self.__str__()


class Combinations:
    """
    Builds combinations of SBML model using antimony

    Create every combination of core hypothesis and extension hypotheses and creates
    SBML models using antimony from the tellurium package.

    :py:class:`Combinations` is designed to be subclassed. The necessary user input
    is given by overriding core functions and providing hypothesis extensions.

    The following methods must be implemented:

        * :py:meth:`core__reactions`
        * :py:meth:`core__parameters`
        * :py:meth:`core__variables`

    However the following methods are optional:

        * :py:meth:`core__functions`
        * :py:meth:`core__events`
        * :py:meth:`core__units`

    Each of these methods should return a valid antimony string since these strings
    are used to build up a full antimony model.

    Extension hypotheses are added by adding methods to your subclass that
    begin with `extension_hypothesis__`. Any method that begins with `extension_hypothesis__`
    will be picked up and used to combinatorially build sbml models.

    Any `extension_hypothesis__` method should return an instance of the
    :py:class:`HypothesisExtension` class, which is merely a container for
    some needed information.

    Extension Hypotheses can operate in either `additive` or `replace` mode,
    depending on how the models should be combined. `additive` is simpler. An extension
    hypothesis is additive when your reaction doesn't override another, or make
    another reaction superflous. Examples of such instances might be
    when adding a mass action reaction to a preexisting set of mass action
    reactions.

    `replace` mode on the other hand should be used when your
    reaction should be used *instead* of another reaction.

    Examples
    --------

    class MyCombModel(Combinations):

        # no __init__ is necessary as we use the __init__ from parent class

        def core__functions(self):
            return ''' '''

        def core__variables(self):
            return '''
compartment Cell;
var A in Cell;
var pA in Cell;
var B in Cell;
var pB in Cell;
var C in Cell;
var pC in Cell;

const S in Cell
'''

        def core__reactions(self):
            return '''
R1f: A -> pA; k1f*A*S;
R2f: B -> pB; k2f*B*A;
R3f: C -> pC; k3f*C*B;
'''

        def core__parameters(self):
            return '''
k1f    = 0.1;
k2f    = 0.1;
k3f    = 0.1;

k2b    = 0.1;
k3b    = 0.1;
VmaxB  = 0.1;
kmB    = 0.1;
VmaxA  = 0.1;
kmA    = 0.1;
k4     = 0.1;
'''

        def core__units(self):
            return None  # Not needed for now

        def core__events(self):
            return None  # No events needed

        def extension_hypothesis__additive1(self):
            return HypothesisExtension(
                name='AdditiveReaction1',
                reaction='pB -> B',
                rate_law='k2b * pB',
                mode='additive',
                to_replace=None,  # not needed for additive mode
            )

        def extension_hypothesis__additive2(self):
            return HypothesisExtension(
                name='AdditiveReaction2',
                reaction='pC -> C',
                rate_law='k3b * C',
                mode='additive',
                to_replace=None,  # not needed for additive mode
            )

        def extension_hypothesis__replace_reaction(self):
            return HypothesisExtension(
                name='ReplaceReaction',
                reaction='pB -> B',
                rate_law='VmaxB * pB / (kmB + pB)',
                mode='replace',
                to_replace='R2f',  # name of reaction we want to replace
            )

        def extension_hypothesis__feedback1(self):
            return HypothesisExtension(
                name='Feedback1',
                reaction='pA -> A',
                rate_law='VmaxA * pA / (kmA + pA)',
                mode='additive',
                to_replace=None,  # name of reaction we want to replace
            )

        def extension_hypothesis__feedback2(self):
            return HypothesisExtension(
                name='Feedback2',
                reaction='pA -> A',
                rate_law='k4 * pA',  # mass action variant
                mode='additive',
                to_replace=None,  # name of reaction we want to replace
            )

    Now that we have built a Combinations subclass we
    can use it as follows:

    .. code-block:
        print(self.c)
        print(len(self.c))
        print(self.c.to_list())

    Produces:

        MyCombModel(topology=0)
        24
        [MyCombModel(topology=0), MyCombModel(topology=1), MyCombModel(topology=2), MyCombModel(topology=3), MyCombModel(topology=4), MyCombModel(topology=5), MyCombModel(topology=6), MyCombModel(topology=7), MyCombModel(topology=8), MyCombModel(topology=9), MyCombModel(topology=10), MyCombModel(topology=11), MyCombModel(topology=12), MyCombModel(topology=13), MyCombModel(topology=14), MyCombModel(topology=15), MyCombModel(topology=16), MyCombModel(topology=17), MyCombModel(topology=18), MyCombModel(topology=19), MyCombModel(topology=20), MyCombModel(topology=21), MyCombModel(topology=22), MyCombModel(topology=23)]


        for i in self.c:
            print(i)

        for i, model in self.c.items():
            print(i, model)

        first_model = self.c[0]

        print(first_model)
        print(first_model.to_antimony())

        rr = first_model.to_tellurium()
        print(rr.simulate(0, 10, 11))








    """

    def __init__(self, directory, mutually_exclusive_reactions=[]):
        self.mutually_exclusive_reactions = mutually_exclusive_reactions
        if self.mutually_exclusive_reactions is not None:
            if not isinstance(self.mutually_exclusive_reactions, list):
                raise TypeError('expecting list but got {}'.format(type(self.mutually_exclusive_reactions)))
            for i in self.mutually_exclusive_reactions:
                if not isinstance(i, tuple):
                    raise TypeError('expecting tuple but got {}'.format(type(self.mutually_exclusive_reactions)))


        self._topology = 0
        self.problem_directory = directory
        if not os.path.isdir(self.problem_directory):
            os.makedirs(self.problem_directory)

        self.cps_file = os.path.join(self.topology_dir, 'Topology{}'.format(self.topology))

        # dict of reactions that vary with topologies and another dict with corresponding hypothesis names
        self.model_variant_reactions, self.topology_names = self._model_variant_reactions()

        # self.model_specific_reactions = self._assembel_model_reactions()[self.topology]

    def __str__(self):
        return "{}(topology={})".format(self.__class__.__name__, self.topology)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        """
        Subtract 1 for 0 indexed python
        :return:
        """
        return len(list(self._get_combinations()))

    def __iter__(self):
        return self

    def __next__(self):

        if self.topology < len(self):
            top = self[self.topology]
            self.topology += 1
            return top
        else:
            self.topology = 0           # turn back to 0 for looping again
            raise StopIteration

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError('"item" should be of type int. Got "{}" instead'.format(type(item)))

        self.topology = item
        return self

    def to_list(self):
        return [deepcopy(self[i]) for i in range(len(self))]

    def items(self):
        return [(i, deepcopy(self[i])) for i in range(len(self))]

    def _model_variant_reactions(self):
        """
        Get all methods that begin with 'extension_hypothesis' and return their values in a dict[number] = reaction_string

        This assembles the reactions that are not in every model and will later be combinatorially combined with the
        core model.

        Returns:

        """
        hypothesis_reactions = []
        hypothesis_reaction_names = []
        for i in dir(self):
            if i.startswith('extension_hypothesis__'):
                hypothesis_reactions.append(getattr(self, i)())
                hypothesis_reaction_names.append(i.replace('extension_hypothesis__', ''))

        dct = OrderedDict()
        names = OrderedDict()
        for i in range(len(hypothesis_reactions)):
            dct[i] = hypothesis_reactions[i]
            names[i] = hypothesis_reaction_names[i]
        return dct, names

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, new):
        assert isinstance(new, int)
        self._topology = new

    @property
    def model_selection_dir(self):
        d = os.path.join(self.problem_directory, 'ModelSelection')
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def topology_dir(self):
        d = os.path.join(self.model_selection_dir, 'Topology{}'.format(self.topology))
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def fit_dir(self):
        d = os.path.join(self.topology_dir, 'Fit{}'.format(self.fit))
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def graphs_dir(self):
        d = os.path.join(self.fit_dir, 'Graphs')
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def time_course_graphs(self):
        d = os.path.join(self.graphs_dir, 'TimeCourseSimulations')
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def copasi_file(self):
        return os.path.join(self.fit_dir, 'topology{}.cps'.format(self.topology))

    def list_topologies(self):
        topologies = OrderedDict()
        comb = self._get_combinations()

        for i in comb:
            if i == ():
                topologies[i] = 'Null'
            else:
                topologies[i] = '__'.join([self.topology_names[x].strip() for x in i])
        # print(topologies)
        df = pandas.DataFrame(topologies, index=['Topology']).transpose().reset_index(drop=True)
        df.index.name = 'ModelID'
        return df

    def to_tellurium(self):
        return te.loada(self.to_antimony())

    def to_antimony(self):
        return self._build_antimony()

    def get_all_parameters_as_list(self):
        all_parameters = self.core__parameters().split('\n')
        all_parameters = [i.strip() for i in all_parameters]
        all_parameters = [re.findall('^\w+', i) for i in all_parameters]
        all_parameters = [i for i in all_parameters if i != []]
        all_parameters = [i[0] for i in all_parameters]
        return all_parameters

    def get_hypotheses(self):
        return self.list_topologies().loc[self.topology][0].split('__')

    def get_reaction_names(self):
        reactions = self.core__reactions().split('\n')
        reactions = [i.strip() for i in reactions]
        reactions = [i for i in reactions if i]
        reactions = [i for i in reactions if not i.startswith('//')]
        names = [re.findall('(.*):', i)[0] for i in reactions]
        return names

    def _get_combinations(self):
        # convert mutually exclusive reactions to numerical value
        mut_excl_list = []
        for mi1, mi2 in self.mutually_exclusive_reactions:
            l2 = []
            for k, v in self.model_variant_reactions.items():
                mi1_match = re.findall('^'+mi1, str(v))
                mi2_match = re.findall('^'+mi2, str(v))

                if mi1_match:
                    l2.append(k)

                if mi2_match:
                    l2.append(k)

            if len(l2) == 1:
                raise ValueError('Cannot have a single reaction '
                                 'in a mutually exclusive pair. Please'
                                 'check that all reactions mentioned '
                                 'in the `mutually_exclusive_reactions` argument'
                                 'actually exist. ')
            mut_excl_list.append(l2)

        # Now that we have the list of MI reactions corresponding to integers, tackle finding the subset
        perm_list = []
        for i in range(len(self.model_variant_reactions)):
            perm_list += [j for j in combinations(range(len(self.model_variant_reactions)), i)]
        # we now need to remove any reaction that contains both parts of a mutually exclusive reaction couple
        perm_list2 = []
        for model_comb in perm_list:
            # when we have no mutually exclusive reactions
            if not mut_excl_list:
                perm_list2.append(model_comb)
            else:
                # when we have mutually exclusive reactions, filter them out
                for mi1, mi2 in mut_excl_list:
                    if mi1 in model_comb and mi2 in model_comb:
                        continue
                    perm_list2.append(model_comb)
        return perm_list2

    def _build_reactions(self):
        """
        Build reactions using two mechanisms. 1) additive. When a HypothesisExtension class is marked as
        additive we can simply add the reaction to the bottom of the list of reactions. 2) replace. Alternatively
        we can replace an existing reaction with the hypothesis
        Returns:

        """
        reactions = self.core__reactions().split('\n')
        reactions = [i.strip() for i in reactions]
        # print(reactions)
        # get additional reactions for current topology

        hypotheses_needed = self._get_combinations()[self._topology]
        hypotheses_needed = [self.model_variant_reactions[i] for i in hypotheses_needed]
        replacements = [i.to_replace for i in hypotheses_needed]
        s = ''
        for reaction in reactions:
            ## reaction name is always the first word, without the colon
            reaction_name = re.findall('^\w+', reaction)

            if reaction_name == []:
                s += '\t\t' + reaction + '\n'
                # continue

            elif reaction_name[0] in replacements:
                # get index of the reaction we want to replace
                idx = replacements.index(reaction_name[0])
                replacement_reaction = hypotheses_needed[idx]
                s += '\t\t' + str(replacement_reaction) + '\n'
            elif reaction_name[0] not in replacements:
                s += '\t\t' + reaction + '\n'
            else:
                raise ValueError('This should not happen')

        # now add the additional extention hypotheses marked as additive
        for i in hypotheses_needed:
            if i.mode == 'additive':
                s += str(i) + '\n'
        return s

    def _build_antimony(self, best_parameters=False):
        """

        :param best_parameters: If False, use default parameters. If
            True, use the best parameters from current fit dir. If a string,
            then it is a parameter set as antimony string
        :return:
        """
        s = ''
        s += self.core__functions()
        s += 'model {}Topology{}'.format(self.__class__.__name__, self.topology)
        s += self.core__variables()
        s += self._build_reactions()

        if best_parameters is False:
            s += self.core__parameters()
        elif best_parameters is True:
            s += self.get_best_model_parameters_as_antimony()

        else:
            raise ValueError
        if self.core__events():
           s += self.core__events()
        if self.core__units():
            s += self.core__units()
        s += "\nend"

        # we now need to remove any global parameters that are not used in the current model topology
        #todo find a better solution for this bit
        exclude_list = ['Cell']  # we want to keep these

        for useless_parameter in self.get_all_parameters_as_list():
            if useless_parameter not in self._build_reactions():
                if useless_parameter not in exclude_list:
                    s = re.sub(useless_parameter + '.*\n', '', s)
        return s

    def _default_parameter_set_as_dict(self):
        string = self.core__parameters()
        strings = string.split('\n')
        dct = OrderedDict()
        for s in strings:
            if s.strip() == '':
                continue
            if ':=' in s:
                k, v = s.split(':=')
            elif '=' in s:
                k, v = s.split('=')

            k = k.strip()
            v = v.replace(';', '')
            try:
                dct[k] = float(v)
            except ValueError:
                dct[k] = v

        return dct

    def core__functions(self):
        return None

    def core__variables(self):
        raise NotImplementedError("You must define your constants, variables and their compartments "
                                  "by defining a `core__variables` method")

    def core__reactions(self):
        raise NotImplementedError('You must define a core set of reactions using the `core__reactions` '
                                  'method.')

    def core__parameters(self):
        raise NotImplementedError('You must define a parameter set (ICs, kinetic parameters, global '
                                  'quantities, compartment volumes) using the `core__parameters` method')

    def core__events(self):
        """
        D = 0
        T = 1
        AZD at t=1.25 == 2
        AZD at t=24 == 3
        AZD at t=48 == 4
        AZD at t=72 == 5
        MK2206 at t=1.25 == 6
        MK2206 at t=24 == 7
        MK2206 at t=48 == 8
        MK2206 at t=72 == 9
        MK2206 and AZD at t=24  == 10
        MK2206 and AZD at t=48 == 11
        MK2206 and AZD at t=72 == 12

        :return:
        """
        return None

    def core__units(self):
        return None



