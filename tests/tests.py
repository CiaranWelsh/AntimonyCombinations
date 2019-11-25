import unittest

from antimony_combinations.antimony_combinations import Combinations, HypothesisExtension
import os
import glob
from shutil import rmtree


class TearDown(unittest.TestCase):

    def tearDown(self) -> None:
        test_dir = os.path.dirname(__file__)
        topology_dirs = glob.glob(os.path.join(test_dir, 'Topology*'))
        for top in topology_dirs:
            if os.path.isdir(top):
                rmtree(top)


class CombinationsTestNoMutualExclusivitity(TearDown):
    class TestCombinationModel(Combinations):

        def core__functions(self):
            return """
                function MM(km, Vmax, S)
                        Vmax * S / (km + S)
                    end

                    function MMWithKcat(km, kcat, S, E)
                        kcat * E * S / (km + S)
                    end


                    function NonCompetitiveInhibition(km, ki, Vmax, n, I, S)
                        Vmax * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end

                    function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
                        kcat * E * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end

                    function NonCompetitiveInhibitionWithKcatAndExtraActivator(km, ki, kcat, E1, E2, n, I, S)
                        kcat * E1 * E2 * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end


                    function MA1(k, S)
                        k * S
                    end

                    function MA2(k, S1, S2)
                        k * S1 * S2
                    end

                    function MA1Mod(k, S, M)
                        k * S * M
                    end

                    function MA2Mod(k, S1, S2, M)
                        k * S1 * S2 * M
                    end

                    function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
                        kcat * E * S / (km + S + ((km * I )/ ki)  )
                    end    

                    function CompetitiveInhibition(Vmax, km, ki, I, S)
                        Vmax * S / (km + S + ((km * I )/ ki)  )
                    end

                    function Hill(km, kcat, L, S, h)
                        kcat * L * (S / km)^h  /   1 + (S / km)^h 
                    end
                """

        def core__variables(self):
            """

            :return:
            """
            return """
                compartment Cell = 1.0

                var Smad2           in Cell  
                var pSmad2          in Cell  
                var Erk             in Cell
                var pErk            in Cell  
                var Akt             in Cell
                var pAkt            in Cell  
                var S6K             in Cell
                var pS6K            in Cell  

                const TGFb             in Cell
                const AZD              in Cell
                const GrowthFactors    in Cell
                const MK2206           in Cell
                const Everolimus       in Cell"""

        def core__reactions(self):
            return """
                //TGFb module
                TGFbR1: Smad2 => pSmad2 ; _kSmad2PhosByTGFb*Smad2*TGFb;
                TGFbR2: pSmad2 => Smad2 ; _kSmad2Dephos*pSmad2;

                //MAPK module
                MAPKR1: Erk => pErk ; kErkPhosByGF*Erk*GrowthFactors;
                MAPKR2: Erk => pErk ; CompetitiveInhibitionWithKcat(_kErkPhosByTGFb_km, _kErkPhosByTGFb_ki, _kErkPhosByTGFb_kcat, TGFb, AZD, Erk);     //(km, ki, kcat, E, I, S)
                MAPKR3: pErk => Erk ; _kErkDephos*pErk;

                //Akt Module
                PI3KR1: Akt => pAkt ; kAktPhosByGF*Akt*GrowthFactors; 
                PI3KR2: Akt => pAkt ; NonCompetitiveInhibitionWithKcat(_kAktPhosByTGFb_km, _kAktPhosByTGFb_km, _kAktPhosByTGFb_kcat, TGFb, 1, MK2206, Akt);  //(km, ki, kcat, E, n, I, S)
                PI3KR3: pAkt => Akt  ; _kAktDephos*pAkt*pS6K;
                PI3KR4: S6K => pS6K ; CompetitiveInhibitionWithKcat(_kS6KPhosByAkt_km, _kS6KPhosByAkt_ki, _kS6KPhosByAkt_kcat, pAkt, Everolimus, S6K); //(km, ki, kcat, E, I, S)
                PI3KR5: pS6K => S6K ; _kS6KDephos*pS6K;

                // Cross talk reactions
            """

        def core__parameters(self):
            return """        

                Akt = 2836.497890686395;
                Erk = 4.720382773676701;
                S6K = 1.5539786683703546e-05;
                Smad2 = 574.6540897866457;
                pAkt = 249.40506182880972;
                pErk = 1.1875677016678898;
                pS6K = 3.943771884555902;
                pSmad2 = 0.0005925125167481985;
                Cell = 1.0;
                AZD = 0.0;
                Everolimus = 0.0;
                ExperimentIndicator = 0.0;
                GrowthFactors = 1.0;
                MK2206 = 0.0;
                TGFb = 0.005;
                _ErkActivateS6K = 0.5146853452098179;
                _kAktActivateErk = 0.00779857827861191;
                _kAktDephos = 2.817381812765289;
                _kAktPhosByTGFb_kcat = 1.9321586563819613;
                _kAktPhosByTGFb_km = 1e-06;
                _kAktPhosSmad2_kcat = 250.09307129734577;
                _kAktPhosSmad2_ki = 2.3711297696305094;
                _kAktPhosSmad2_km = 159849.61877127536;
                _kErkDephos = 0.8120002369435844;
                _kErkPhosByTGFb_kcat = 7349.566965163551;
                _kErkPhosByTGFb_ki = 999944.5901705284;
                _kErkPhosByTGFb_km = 83900.71885922311;
                _kS6KDephos = 1.6547141713211548;
                _kS6KPhosByAkt_kcat = 2393.799210262596;
                _kS6KPhosByAkt_ki = 1.000002108620802e-06;
                _kS6KPhosByAkt_km = 181.24417678017062;
                _kSmad2Dephos = 0.18123319271852234;
                kAktPhosByGF = 0.1;
                kErkPhosByGF = 0.1;

        		"""

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
            return """
                // events in all simulations
                SerumStarveRemoveTGFb: at (time>70.25): TGFb=0.00005;
                SerumStarveRemoveGrowthFactors: at (time>70.25): GrowthFactors=0.005;

                // these events are dependent on the experiment indicated by the ExperimentIndicator Variable
                AddTGFb:        at (time>71.25  and ExperimentIndicator >  0):   TGFb=1;
                AddAZD_1_25:    at (time>70.75  and ExperimentIndicator == 2):   AZD=1;
                AddAZD_24:      at  (time>48    and ExperimentIndicator == 3):   AZD=1;
                AddAZD_48:      at  (time>24    and ExperimentIndicator == 4):   AZD=1;
                AddAZD_72:      at  (time>0     and ExperimentIndicator == 5):   AZD=1;
                AddMK_1_25:     at (time>70.75  and ExperimentIndicator == 6):   MK2206=1;
                AddMK_24:       at (time>48     and ExperimentIndicator == 7):   MK2206=1;
                AddMK_48:       at (time>24     and ExperimentIndicator == 8):   MK2206=1;
                AddMK_72:       at (time>0      and ExperimentIndicator == 9):   MK2206=1;
                AddAZDAndMK_24: at (time>48     and ExperimentIndicator == 10):  MK2206=1, AZD=1;
                AddAZDAndMK_48: at (time>24     and ExperimentIndicator == 11):  MK2206=1, AZD=1;
                AddAZDAndMK_72: at (time>0      and ExperimentIndicator == 12):  MK2206=1, AZD=1;
                """

        def core__units(self):
            return """
                unit volume = 1 litre;
                unit time_unit = 3600 second;
                unit substance = 1e-9 mole;
                """

        def extension_hypothesis__AktActivateSmad2ErkInhibit(self):
            """
            This reaction must replace:
                TGFbR1: Smad2 => pSmad2 ; _kSmad2PhosByTGFb*Smad2*TGFb;
            Args:
                type:
                replacement_reaction:

            Returns:

            """
            return HypothesisExtension(
                name='CrossTalkR1',
                reaction='Smad2 => pSmad2',
                rate_law='NonCompetitiveInhibitionWithKcatAndExtraActivator(_kAktPhosSmad2_km, _kAktPhosSmad2_ki, _kAktPhosSmad2_kcat, TGFb, pAkt, 1, pErk, Smad2)',
                mode='replace',
                to_replace='TGFbR1'
            )

        def extension_hypothesis__ErkActivateSmad2AktInhibit(self):
            return HypothesisExtension(
                name='CrossTalkR2',
                reaction='Smad2 => pSmad2',
                rate_law='NonCompetitiveInhibitionWithKcatAndExtraActivator(_kErkPhosSmad2_km, _kErkPhosSmad2_ki, _kErkPhosSmad2_kcat, TGFb, pErk, 1, pAkt, Smad2);  //(km, ki, kcat, E, n, I, S)',
                mode='replace',
                to_replace='TGFbR1'
            )

        def extension_hypothesis__pAktActivateErk(self):
            return HypothesisExtension(
                name='CrossTalkR4',
                reaction='Erk => pErk',
                rate_law='_kAktActivateErk*Erk*pAkt',
                mode='additive',
                to_replace=None
            )

        def extension_hypothesis__S6KActivateErk(self):
            return HypothesisExtension(
                name='CrossTalkR5',
                reaction='Erk => pErk',
                rate_law='_kS6KActivateErk*Erk*pS6K',
                mode='additive',
                to_replace=None
            )

        def extension_hypothesis__ErkActivatesS6K(self):
            return HypothesisExtension(
                name='CrossTalkR6',
                reaction='S6K => pS6K',
                rate_law='_ErkActivateS6K*pErk*S6K',
                mode='additive',
                to_replace=None
            )

    def setUp(self) -> None:
        self.dire = os.path.dirname(__file__)
        self.c = self.TestCombinationModel(directory=self.dire)

    def test_len(self):
        expected = 31
        actual = len(self.c)
        self.assertEqual(expected, actual)

    def test_getitem(self):
        mod4 = self.c[4]
        self.assertEqual(mod4.topology, 4)

    def test_get_hypotheses(self):
        mod4 = self.c[16]
        expected = ['AktActivateSmad2ErkInhibit', 'ErkActivateSmad2AktInhibit', 'ErkActivatesS6K']
        actual = mod4.get_hypotheses()
        self.assertEqual(expected, actual)

    def test_iterater(self):
        count = 0
        for id in self.c:
            count += 1
        self.assertEqual(count, len(self.c))

    def test_topology_counter_put_back_to_0(self):
        count = 0
        for id in self.c:
            count += 1

        self.assertEqual(self.c.topology, 0)

    def test_items_method(self):
        l = self.c.items()
        actual = sorted([i[0] for i in self.c.items()])
        expected = sorted(list(range(31)))
        self.assertEqual(expected, actual)

    def test_to_list(self):
        l = self.c.to_list()
        actual = [i.topology for i in l]
        expected = list(range(31))
        self.assertEqual(expected, actual)

    def test_mutually_exclusive_models(self):
        pass

    def test_list_constructor(self):
        print(list(self.c))

    def test_get_reactioin_names(self):
        expected = ['TGFbR1', 'TGFbR2', 'MAPKR1', 'MAPKR2', 'MAPKR3', 'PI3KR1', 'PI3KR2', 'PI3KR3', 'PI3KR4', 'PI3KR5']
        actual = self.c.get_reaction_names()
        self.assertEqual(expected, actual)

    def test_slice(self):
        expected = [0, 2, 4, 6, 8]
        actual = [i.topology for i in self.c[:10:2]]
        self.assertEqual(expected, actual)

    def test_subset_by_list(self):
        actual = [4, 9]
        expected = [i.topology for i in self.c[actual]]
        self.assertEqual(expected, actual)

    def test_subset_by_tuple(self):
        actual = (7, 5)
        expected = [i.topology for i in self.c[actual]]
        self.assertEqual(expected, list(actual))


class CombinationsTestWithMutualExclusivitity(TearDown):
    class TestCombinationModel(Combinations):

        def core__functions(self):
            return """
                function MM(km, Vmax, S)
                        Vmax * S / (km + S)
                    end

                    function MMWithKcat(km, kcat, S, E)
                        kcat * E * S / (km + S)
                    end


                    function NonCompetitiveInhibition(km, ki, Vmax, n, I, S)
                        Vmax * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end

                    function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
                        kcat * E * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end

                    function NonCompetitiveInhibitionWithKcatAndExtraActivator(km, ki, kcat, E1, E2, n, I, S)
                        kcat * E1 * E2 * S / ( (km + S) * (1 + (I / ki)^n ) )
                    end


                    function MA1(k, S)
                        k * S
                    end

                    function MA2(k, S1, S2)
                        k * S1 * S2
                    end

                    function MA1Mod(k, S, M)
                        k * S * M
                    end

                    function MA2Mod(k, S1, S2, M)
                        k * S1 * S2 * M
                    end

                    function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
                        kcat * E * S / (km + S + ((km * I )/ ki)  )
                    end    

                    function CompetitiveInhibition(Vmax, km, ki, I, S)
                        Vmax * S / (km + S + ((km * I )/ ki)  )
                    end

                    function Hill(km, kcat, L, S, h)
                        kcat * L * (S / km)^h  /   1 + (S / km)^h 
                    end
                """

        def core__variables(self):
            """

            :return:
            """
            return """
                compartment Cell = 1.0

                var Smad2           in Cell  
                var pSmad2          in Cell  
                var Erk             in Cell
                var pErk            in Cell  
                var Akt             in Cell
                var pAkt            in Cell  
                var S6K             in Cell
                var pS6K            in Cell  

                const TGFb             in Cell
                const AZD              in Cell
                const GrowthFactors    in Cell
                const MK2206           in Cell
                const Everolimus       in Cell"""

        def core__reactions(self):
            return """
                //TGFb module
                TGFbR1: Smad2 => pSmad2 ; _kSmad2PhosByTGFb*Smad2*TGFb;
                TGFbR2: pSmad2 => Smad2 ; _kSmad2Dephos*pSmad2;

                //MAPK module
                MAPKR1: Erk => pErk ; kErkPhosByGF*Erk*GrowthFactors;
                MAPKR2: Erk => pErk ; CompetitiveInhibitionWithKcat(_kErkPhosByTGFb_km, _kErkPhosByTGFb_ki, _kErkPhosByTGFb_kcat, TGFb, AZD, Erk);     //(km, ki, kcat, E, I, S)
                MAPKR3: pErk => Erk ; _kErkDephos*pErk;

                //Akt Module
                PI3KR1: Akt => pAkt ; kAktPhosByGF*Akt*GrowthFactors; 
                PI3KR2: Akt => pAkt ; NonCompetitiveInhibitionWithKcat(_kAktPhosByTGFb_km, _kAktPhosByTGFb_km, _kAktPhosByTGFb_kcat, TGFb, 1, MK2206, Akt);  //(km, ki, kcat, E, n, I, S)
                PI3KR3: pAkt => Akt  ; _kAktDephos*pAkt*pS6K;
                PI3KR4: S6K => pS6K ; CompetitiveInhibitionWithKcat(_kS6KPhosByAkt_km, _kS6KPhosByAkt_ki, _kS6KPhosByAkt_kcat, pAkt, Everolimus, S6K); //(km, ki, kcat, E, I, S)
                PI3KR5: pS6K => S6K ; _kS6KDephos*pS6K;

                // Cross talk reactions
            """

        def core__parameters(self):
            return """        

                Akt = 2836.497890686395;
                Erk = 4.720382773676701;
                S6K = 1.5539786683703546e-05;
                Smad2 = 574.6540897866457;
                pAkt = 249.40506182880972;
                pErk = 1.1875677016678898;
                pS6K = 3.943771884555902;
                pSmad2 = 0.0005925125167481985;
                Cell = 1.0;
                AZD = 0.0;
                Everolimus = 0.0;
                ExperimentIndicator = 0.0;
                GrowthFactors = 1.0;
                MK2206 = 0.0;
                TGFb = 0.005;
                _ErkActivateS6K = 0.5146853452098179;
                _kAktActivateErk = 0.00779857827861191;
                _kAktDephos = 2.817381812765289;
                _kAktPhosByTGFb_kcat = 1.9321586563819613;
                _kAktPhosByTGFb_km = 1e-06;
                _kAktPhosSmad2_kcat = 250.09307129734577;
                _kAktPhosSmad2_ki = 2.3711297696305094;
                _kAktPhosSmad2_km = 159849.61877127536;
                _kErkDephos = 0.8120002369435844;
                _kErkPhosByTGFb_kcat = 7349.566965163551;
                _kErkPhosByTGFb_ki = 999944.5901705284;
                _kErkPhosByTGFb_km = 83900.71885922311;
                _kS6KDephos = 1.6547141713211548;
                _kS6KPhosByAkt_kcat = 2393.799210262596;
                _kS6KPhosByAkt_ki = 1.000002108620802e-06;
                _kS6KPhosByAkt_km = 181.24417678017062;
                _kSmad2Dephos = 0.18123319271852234;
                kAktPhosByGF = 0.1;
                kErkPhosByGF = 0.1;

        		"""

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
            return """
                // events in all simulations
                SerumStarveRemoveTGFb: at (time>70.25): TGFb=0.00005;
                SerumStarveRemoveGrowthFactors: at (time>70.25): GrowthFactors=0.005;

                // these events are dependent on the experiment indicated by the ExperimentIndicator Variable
                AddTGFb:        at (time>71.25  and ExperimentIndicator >  0):   TGFb=1;
                AddAZD_1_25:    at (time>70.75  and ExperimentIndicator == 2):   AZD=1;
                AddAZD_24:      at  (time>48    and ExperimentIndicator == 3):   AZD=1;
                AddAZD_48:      at  (time>24    and ExperimentIndicator == 4):   AZD=1;
                AddAZD_72:      at  (time>0     and ExperimentIndicator == 5):   AZD=1;
                AddMK_1_25:     at (time>70.75  and ExperimentIndicator == 6):   MK2206=1;
                AddMK_24:       at (time>48     and ExperimentIndicator == 7):   MK2206=1;
                AddMK_48:       at (time>24     and ExperimentIndicator == 8):   MK2206=1;
                AddMK_72:       at (time>0      and ExperimentIndicator == 9):   MK2206=1;
                AddAZDAndMK_24: at (time>48     and ExperimentIndicator == 10):  MK2206=1, AZD=1;
                AddAZDAndMK_48: at (time>24     and ExperimentIndicator == 11):  MK2206=1, AZD=1;
                AddAZDAndMK_72: at (time>0      and ExperimentIndicator == 12):  MK2206=1, AZD=1;
                """

        def core__units(self):
            return """
                unit volume = 1 litre;
                unit time_unit = 3600 second;
                unit substance = 1e-9 mole;
                """

        def extension_hypothesis__AktActivateSmad2ErkInhibit(self):
            """
            This reaction must replace:
                TGFbR1: Smad2 => pSmad2 ; _kSmad2PhosByTGFb*Smad2*TGFb;
            Args:
                type:
                replacement_reaction:

            Returns:

            """
            return HypothesisExtension(
                name='CrossTalkR1',
                reaction='Smad2 => pSmad2',
                rate_law='NonCompetitiveInhibitionWithKcatAndExtraActivator(_kAktPhosSmad2_km, _kAktPhosSmad2_ki, _kAktPhosSmad2_kcat, TGFb, pAkt, 1, pErk, Smad2)',
                mode='replace',
                to_replace='TGFbR1'
            )

        def extension_hypothesis__ErkActivateSmad2AktInhibit(self):
            return HypothesisExtension(
                name='CrossTalkR2',
                reaction='Smad2 => pSmad2',
                rate_law='NonCompetitiveInhibitionWithKcatAndExtraActivator(_kErkPhosSmad2_km, _kErkPhosSmad2_ki, _kErkPhosSmad2_kcat, TGFb, pErk, 1, pAkt, Smad2);  //(km, ki, kcat, E, n, I, S)',
                mode='replace',
                to_replace='TGFbR1'
            )

        def extension_hypothesis__pAktActivateErk(self):
            return HypothesisExtension(
                name='CrossTalkR4',
                reaction='Erk => pErk',
                rate_law='_kAktActivateErk*Erk*pAkt',
                mode='additive',
                to_replace=None
            )

        def extension_hypothesis__S6KActivateErk(self):
            return HypothesisExtension(
                name='CrossTalkR5',
                reaction='Erk => pErk',
                rate_law='_kS6KActivateErk*Erk*pS6K',
                mode='additive',
                to_replace=None
            )

        def extension_hypothesis__ErkActivatesS6K(self):
            return HypothesisExtension(
                name='CrossTalkR6',
                reaction='S6K => pS6K',
                rate_law='_ErkActivateS6K*pErk*S6K',
                mode='additive',
                to_replace=None
            )

    def setUp(self) -> None:
        self.dire = os.path.dirname(__file__)
        self.c = self.TestCombinationModel(
            mutually_exclusive_reactions=[
                ('CrossTalkR1', 'CrossTalkR2'),
            ], directory=self.dire
        )

    def test_get_combinations(self):
        self.c._get_combinations()

    def test_len(self):
        expected = 24
        actual = len(self.c)
        self.assertEqual(expected, actual)

    def test_getitem(self):
        mod4 = self.c[4]
        self.assertEqual(mod4.topology, 4)


class AnotherExampleTests(TearDown):
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

    def setUp(self) -> None:
        directory = os.path.dirname(__file__)
        self.c = self.MyCombModel(mutually_exclusive_reactions=[
            ('Feedback1', 'Feedback2')
        ], directory=directory)

    def test__output_used_in_docs(self):
        """
        Keep for now. You may want to update the docs.
        Returns:

        """
        print(self.c)
        print(len(self.c))
        print(self.c.to_list()[:4])

        for i in self.c[:3]:
            print(i)

        for i, model in self.c.items()[:3]:
            print(i, model)

        first_model = self.c[0]

        print(first_model)
        print(first_model.to_antimony())

        rr = first_model.to_roadrunner()
        print(rr)
        print(rr.simulate(0, 10, 11))

        print(self.c.get_topologies())

        print(self.c.to_copasi())


if __name__ == '__main__':
    unittest.main()
