import site, os, glob
import pandas, numpy
import re
import tellurium as te
import roadrunner as rr
import typing

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
from pycotools3 import tasks, model

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

tuple_list = typing.List[typing.Tuple[typing.AnyStr]]

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

    The following methods must be implemented (see below for an example):

        * :py:meth:`core__reactions`
        * :py:meth:`core__parameters`
        * :py:meth:`core__variables`

    However the following methods are optional:

        * :py:meth:`core__functions`
        * :py:meth:`core__events`
        * :py:meth:`core__units`

    Each of these methods should return a
    `valid antimony string <https://tellurium.readthedocs.io/en/latest/antimony.html>`_,
    since these strings are used to build up a full antimony model.

    Extension hypotheses are added by adding methods to your subclass that
    begin with `extension_hypothesis__`. Any method that begins with `extension_hypothesis__`
    will be picked up and used to combinatorially build sbml models.

    Any `extension_hypothesis__` method should return an instance of the
    :py:class:`HypothesisExtension` class, which is merely a container for
    some needed information.

    .. note::

        Notice the double underscore after `extension_hypothesis`

    Extension Hypotheses can operate in either `additive` or `replace` mode,
    depending on how the models should be combined. `additive` is simpler. An extension
    hypothesis is additive when your reaction doesn't override another, or make
    another reaction superflous. Examples of such instances might be
    when adding a mass action reaction to a preexisting set of mass action
    reactions.

    `replace` mode on the other hand should be used when your
    reaction should be used *instead* of another reaction.


    Examples:

    .. code-block:: python
        :linenos:

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

                S = 1;
                A = 10;
                pA = 0;
                B = 10;
                pB = 0;
                C = 10;
                pC = 0;
                Cell = 1;
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

    >>> project_root = os.path.dirname(__file__)
    >>> c = MyCombModel(mutually_exclusive_reactions=[
    >>>         ('Feedback1', 'Feedback2')
    >>>     ], directory=project_root       # optionally specify project root
    >>> )

    MyCombModel behaves like an iterator, though it doesn't
    store all model topologies on the outset but builds models
    of the fly as the `topology` attribute is incremented. Topology
    always starts on model 0, the core model that
    doesn't have additional hypothesis extensions.

    >>> print(c)
    MyCombModel(topology=0)

    The complete set of model topologies is enumerated by the
    `topology` attribute. The `__len__` method is set to the
    size of this set, accounting for mutually exclusive
    topologies, which is a mechanism for reducing the topology space.

    >>> print(len(c))
    24

    You can pick out any of these topologies using the
    selection operator

    >>> print(c[4])
    MyCombModel(topology=4)

    To see which topologies correspond to which hypothesis
    extensions we can use :py:meth:`antimony_combinations.list_topologies`,
    which returns a pandas.DataFrame.

    >>> c.list_topolgies()
                                                          Topology
    ModelID
    0                                                     Null
    1                                                additive1
    2                                                additive2
    3                                                feedback1
    4                                                feedback2
    5                                         replace_reaction
    6                                     additive1__additive2
    7                                     additive1__feedback1
    8                                     additive1__feedback2
    9                              additive1__replace_reaction
    10                                    additive2__feedback1
    11                                    additive2__feedback2
    12                             additive2__replace_reaction
    13                             feedback1__replace_reaction
    14                             feedback2__replace_reaction
    15                         additive1__additive2__feedback1
    16                         additive1__additive2__feedback2
    17                  additive1__additive2__replace_reaction
    18                  additive1__feedback1__replace_reaction
    19                  additive1__feedback2__replace_reaction
    20                  additive2__feedback1__replace_reaction
    21                  additive2__feedback2__replace_reaction
    22       additive1__additive2__feedback1__replace_reaction
    23       additive1__additive2__feedback2__replace_reaction


    You can extract all topologies into a list using the
    :py:meth:`antimony_combinations.Combinations.to_list`
    method.

    >>> print(c.to_list()[:4])
    [MyCombModel(topology=0),
     MyCombModel(topology=1),
     MyCombModel(topology=2),
     MyCombModel(topology=3)]

    You can iterate over the set of topologies

    >>> for i in c[:3]:
    >>> ... print(i)
    MyCombModel(topology=0)
    MyCombModel(topology=1)
    MyCombModel(topology=2)

    Or use the items method, which is similar to `dict.items()`.

    >>> for i, model in c.items()[:3]:
    >>> ... print(i, model)
    0 MyCombModel(topology=0)
    1 MyCombModel(topology=1)
    2 MyCombModel(topology=2)

    Selecting a single model, we can create an antimony string

    >>> first_model = c[0]
    >>> print(first_model.to_antimony())
    model MyCombModelTopology0
        compartment Cell;
        var A in Cell;
        var pA in Cell;
        var B in Cell;
        var pB in Cell;
        var C in Cell;
        var pC in Cell;
        const S in Cell
        R1f: A -> pA; k1f*A*S;
        R2f: B -> pB; k2f*B*A;
        R3f: C -> pC; k3f*C*B;
        k1f = 0.1;
        k2f = 0.1;
        k3f = 0.1;
        S = 1;
        A = 10;
        pA = 0;
        B = 10;
        pB = 0;
        C = 10;
        pC = 0;
        Cell = 1;
    end

    or a tellurium model

    >>> rr = first_model.to_tellurium()
    >>> print(rr)
    <roadrunner.RoadRunner() {
    'this' : 0x555a52c8cb90
    'modelLoaded' : true
    'modelName' :
    'libSBMLVersion' : LibSBML Version: 5.17.2
    'jacobianStepSize' : 1e-05
    'conservedMoietyAnalysis' : false
    'simulateOptions' :
    < roadrunner.SimulateOptions()
    {
    'this' : 0x555a5309cd00,
    'reset' : 0,
    'structuredResult' : 0,
    'copyResult' : 1,
    'steps' : 50,
    'start' : 0,
    'duration' : 5
    }>,
    'integrator' :
    < roadrunner.Integrator() >
      name: cvode
      settings:
          relative_tolerance: 0.000001
          absolute_tolerance: 0.000000000001
                       stiff: true
           maximum_bdf_order: 5
         maximum_adams_order: 12
           maximum_num_steps: 20000
           maximum_time_step: 0
           minimum_time_step: 0
           initial_time_step: 0
              multiple_steps: false
          variable_step_size: false

    }>

    >>> print(rr.simulate(0, 10, 11))
        time,     [A],     [pA],       [B],    [pB],     [C],    [pC]
     [[    0,      10,        0,        10,       0,      10,       0],
      [    1, 9.04837, 0.951626,   3.86113, 6.13887, 5.27257, 4.72743],
      [    2, 8.18731,  1.81269,   1.63214, 8.36786, 4.07751, 5.92249],
      [    3, 7.40818,  2.59182,  0.748842, 9.25116, 3.64313, 6.35687],
      [    4,  6.7032,   3.2968,  0.370018, 9.62998, 3.45361, 6.54639],
      [    5, 6.06531,  3.93469,  0.195519, 9.80448,  3.3609,  6.6391],
      [    6, 5.48812,  4.51188,  0.109779, 9.89022, 3.31158, 6.68842],
      [    7, 4.96585,  5.03415, 0.0651185, 9.93488,  3.2835,  6.7165],
      [    8, 4.49329,  5.50671, 0.0405951,  9.9594, 3.26657, 6.73343],
      [    9,  4.0657,   5.9343, 0.0264712, 9.97353, 3.25584, 6.74416],
      [   10, 3.67879,  6.32121, 0.0179781, 9.98202, 3.24872, 6.75128]]

    Or an interface to copasi, via `pycotools3 <https://pycotools3.readthedocs.io/en/latest/>`_

    >>> c.to_copasi()
    Model(name=NoName, time_unit=s, volume_unit=l, quantity_unit=mol)

    Which could be used to configure parameter estimations. Currently,
    support for parameter estimation configuration has in COPASI not been included
    but this is planned for the near future.
    """

    def __init__(self,
                 mutually_exclusive_reactions: tuple_list = [],
                 directory: typing.Optional[str] = None) -> None:
        """

        Args:
            mutually_exclusive_reactions:
                An arbitrary length list of tuples of pairs that are names
                of reactions that should never occur together in the same model.
                Defaults to an empty list.
            directory:
                Root directory for analysis. The default is the directory
                containing the script being run or the current working directory
                of the interpreter.
        """

        self.mutually_exclusive_reactions = mutually_exclusive_reactions
        if self.mutually_exclusive_reactions is not None:
            if not isinstance(self.mutually_exclusive_reactions, list):
                raise TypeError('expecting list but got {}'.format(type(self.mutually_exclusive_reactions)))
            for i in self.mutually_exclusive_reactions:
                if not isinstance(i, tuple):
                    raise TypeError('expecting tuple but got {}'.format(type(self.mutually_exclusive_reactions)))

        self._topology = 0
        self.directory = directory

        # for when run from script
        if self.directory is None:
            self.directory = os.path.dirname(__file__)

        # for when run in python interpreter
        if self.directory == '':
            self.directory = os.getcwd()

        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

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
            self.topology = 0  # turn back to 0 for looping again
            raise StopIteration

    def __getitem__(self, item):
        if not isinstance(item, (int, slice,
                                 tuple, list)):
            raise TypeError('"item" should be of type int or slice. Got "{}" instead'.format(type(item)))
        if isinstance(item, int):
            self.topology = item
            return self
        elif isinstance(item, slice):
            required = range(item.start if item.start is not None else 0,
                             item.stop if item.stop is not None else len(self),
                             item.step if item.step is not None else 1)
            return [deepcopy(self[i]) for i in required]
        elif isinstance(item, (tuple, list)):
            for i in item:
                if not isinstance(i, int):
                    raise ValueError('expected an integer for index. Got "{}"'.format(type(i)))
            return [deepcopy(self[i]) for i in item]

    def to_list(self) -> list:
        return [deepcopy(self[i]) for i in range(len(self))]

    def items(self) -> list:
        return [(i, deepcopy(self[i])) for i in range(len(self))]

    def _model_variant_reactions(self) -> typing.Tuple[typing.Dict[int, str], str]:
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
    def topology(self) -> int:
        return self._topology

    @topology.setter
    def topology(self, new) -> None:
        assert isinstance(new, int)
        self._topology = new

    @property
    def topology_dir(self) -> str:
        d = os.path.join(self.directory, 'Topology{}'.format(self.topology))
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def time_course_graphs(self) -> str:
        d = os.path.join(self.topology_dir, 'TimeCourseSimulations')
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    @property
    def copasi_file(self) -> str:
        return os.path.join(self.topology_dir, 'topology{}.cps'.format(self.topology))

    def to_copasi(self) -> model.Model:
        return model.loada(self.to_antimony(), self.copasi_file)

    def list_topologies(self) -> pandas.DataFrame:
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

    def to_tellurium(self) -> rr.ExecutableModel:
        return te.loada(self.to_antimony())

    def to_antimony(self) -> str:
        return self._build_antimony()

    def get_all_parameters_as_list(self) -> typing.List[str]:
        all_parameters = self.core__parameters().split('\n')
        all_parameters = [i.strip() for i in all_parameters]
        all_parameters = [re.findall('^\w+', i) for i in all_parameters]
        all_parameters = [i for i in all_parameters if i != []]
        all_parameters = [i[0] for i in all_parameters]
        return all_parameters

    def get_hypotheses(self) -> typing.List[str]:
        return self.list_topologies().loc[self.topology][0].split('__')

    def get_reaction_names(self) -> typing.List[str]:
        reactions = self.core__reactions().split('\n')
        reactions = [i.strip() for i in reactions]
        reactions = [i for i in reactions if i]
        reactions = [i for i in reactions if not i.startswith('//')]
        names = [re.findall('(.*):', i)[0] for i in reactions]
        return names

    def _get_combinations(self) -> typing.List[typing.Tuple[int]]:
        # convert mutually exclusive reactions to numerical value
        mut_excl_list = []
        for mi1, mi2 in self.mutually_exclusive_reactions:
            l2 = []
            for k, v in self.model_variant_reactions.items():
                mi1_match = re.findall('^' + mi1, str(v))
                mi2_match = re.findall('^' + mi2, str(v))

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

    def _build_reactions(self) -> str:
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
                s += '\t' + reaction + '\n'
                # continue

            elif reaction_name[0] in replacements:
                # get index of the reaction we want to replace
                idx = replacements.index(reaction_name[0])
                replacement_reaction = hypotheses_needed[idx]
                s += '\t' + str(replacement_reaction) + '\n'
            elif reaction_name[0] not in replacements:
                s += '\t' + reaction + '\n'
            else:
                raise ValueError('This should not happen')

        # now add the additional extention hypotheses marked as additive
        for i in hypotheses_needed:
            if i.mode == 'additive':
                s += str(i) + '\n'
        return s

    def _build_antimony(self, best_parameters=False) -> str:
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
        # todo find a better solution for this bit
        exclude_list = ['Cell']  # we want to keep these

        for useless_parameter in self.get_all_parameters_as_list():
            if useless_parameter not in self._build_reactions():
                if useless_parameter not in exclude_list:
                    s = re.sub(useless_parameter + '.*\n', '', s)
        return s

    def _default_parameter_set_as_dict(self) -> typing.Dict[str, float]:
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
        """
        An optional set of functions for use in rate laws. Do not use directly
        but instead override in subclass.

        For example:

        .. code-block::
            :linenos:

            def core__functions(self):
                return '''
                function MichaelisMenten(vmax, km, s)
                    vmax * s / (km + s)
                end
                '''

        Returns (str):

        """
        return None

    def core__variables(self):
        """
        List your variables whilst specifying their compartment.

        Method not to be used directly but overriden in subclass. This is
        a required method.

        Examples:

            .. code-block::
                :linenos:

                def core__variables(self):
                    return '''
                    compartment Cell = 1;
                    var A in Cell;
                    var B in Cell;
                    const S in Cell;
                    '''

        Returns (str):

        """
        raise NotImplementedError("You must define your constants, variables and their compartments "
                                  "by defining a `core__variables` method")

    def core__reactions(self):
        """
        List of core reactions; reactions to be shared among all models.

        Do not use directly as this method is designed to be subclassed. This method
        is required.

        Examples:

            .. code:block::
                :linenos:

                def core__reactions(self):
                    return '''
                    R1: A -> B; k1 * A
                    '''

        Returns (str):

        """
        raise NotImplementedError('You must define a core set of reactions using the `core__reactions` '
                                  'method.')

    def core__parameters(self):
        """
        Parameter list. Do not use directly but over ride in subclass. This method
        is required.

        Examples:

            .. code-block::
                :linenos:

                def core__parameters(self):
                    return '''
                    A = 10;
                    B = 20;

                    k1 = 0.1;
                    hypothesis_extension_parameter1 = 10;   // will be pruned if not in present model.

                    random_global_parameter = 50;
                    '''
        Returns (str):

        """
        raise NotImplementedError('You must define a parameter set (ICs, kinetic parameters, global '
                                  'quantities, compartment volumes) using the `core__parameters` method')

    def core__events(self):
        """
        Antimony events string.

        Do not use directly but override in subclass. Optional method.

        Examples:

            .. code-block::
                :linenos:

                def core__events(self):
                    return '''
                    event1 at t=1.25 == 2
                    '''
        Returns (str):
        """
        return None

    def core__units(self):
        """
        Antimony units string.

        Do not use directly but override in subclass. Optional method.

        Returns:

        """
        return None
