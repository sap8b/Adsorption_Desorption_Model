Option Strict On
Imports System
Imports System.IO

Module MainModule

    Public Class ModelProperties
        Public type_of_simulation As Simulation_Type
        Public number_of_simulation_repetitions, total_time_steps_per_simulation As Integer
        Public temperature As Double
        Public alloy_material As Alloy
        Public reaction_surface As Lattice_Properties
        Public solution As Electrolyte
        Public analytical_model_available As Boolean
        Public reaction_rates_filename As String
        Public coverage_filename As String
        Public charge_transfer_filename As String
        Public final_configuration_filename As String
        Public geometry As Cell_Geometry
        Public computational_box_dimensions As ComputationalCellParameters '() As Integer
        Public reactions() As Kinetics_Types
        Public reaction_list() As Reaction_Process
        Public W_probs() As Double
    End Class

    Sub Main()
        Console.WriteLine("Starting KMC Simulation...")
        Dim KMCModel_Setup As ModelProperties = New ModelProperties With {
            .type_of_simulation = Simulation_Type.Isotherm_Model_Without_Neighbor_Interactions,'Simulation_Type.Isotherm_Model_With_Neighbor_Interaction_Blocking, '
            .number_of_simulation_repetitions = 10, 'Set the number of repetitions of the simulations to run for averaging
            .total_time_steps_per_simulation = 10000, 'Set the length of each simulation in simulation steps.  Time, of course, will vary.
            .temperature = 298,
            .alloy_material = Alloy.Gold,
            .solution = Electrolyte.Water,
            .reactions = New Kinetics_Types() {Kinetics_Types.Adsorption_Basic, Kinetics_Types.Desorption_Basic}, '{Kinetics_Types.Adsorption_WithCharge, Kinetics_Types.Desorption_WithCharge},
            .analytical_model_available = True,
            .reaction_rates_filename = "reaction_rates.dat",
            .coverage_filename = "coverage_with_analytical_solution.dat",
            .charge_transfer_filename = "current_vs_time.dat",
            .final_configuration_filename = "configuration.dat",
            .geometry = Cell_Geometry.Flat 'Cell_Geometry.Rectangular_Parallelpiped
        }
        KMCModel_Setup.reaction_surface = DetermineLatticeProperties(KMCModel_Setup.alloy_material)
        KMCModel_Setup.computational_box_dimensions = Determine_BoxSize(KMCModel_Setup.geometry, KMCModel_Setup.reaction_surface) 'Determine_BoxSize(Cell_Geometry.Flat),'New Integer() {50, 50, 1},

        ReDim KMCModel_Setup.reaction_list(KMCModel_Setup.reactions.Length)
        ReDim KMCModel_Setup.W_probs(KMCModel_Setup.reactions.Length)

        Console.WriteLine("Calculating reaction rates...")
        For r_idx As Integer = 0 To KMCModel_Setup.reactions.Length - 1
            KMCModel_Setup.reaction_list(r_idx) = Determine_Reaction_Rates(KMCModel_Setup.temperature, KMCModel_Setup.reactions(r_idx)) 'New Double() {0.5 * 2, 0.5 * 2}  'This initializes the number of reactions to be considered
            'Console.WriteLine(KMCModel_Setup.reaction_list(r_idx).k_reaction)
        Next
        WriteReactionRates(KMCModel_Setup.reaction_rates_filename, KMCModel_Setup.reaction_list)

        'Console.WriteLine("Calculating transition probabilities...")

        Dim TrackingArrays As TrackingParameters = New TrackingParameters With {
            .plot_sites = New LinkedList(Of ReactionPosition)}
        ReDim TrackingArrays.tau(KMCModel_Setup.number_of_simulation_repetitions, KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.theta(KMCModel_Setup.number_of_simulation_repetitions, KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.charge(KMCModel_Setup.number_of_simulation_repetitions, KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.charge_tau(KMCModel_Setup.number_of_simulation_repetitions, KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.temp_tau(KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.temp_charge(KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.temp_theta(KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.analytical_tau(KMCModel_Setup.total_time_steps_per_simulation)
        ReDim TrackingArrays.analytical_theta(KMCModel_Setup.total_time_steps_per_simulation)

        Console.WriteLine("Beginning repetitions of each simulation...")
        Console.WriteLine()
        For rep_counter As Integer = 0 To KMCModel_Setup.number_of_simulation_repetitions - 1

            TrackingArrays.plot_sites.Clear()

            Select Case KMCModel_Setup.type_of_simulation
                Case Simulation_Type.Isotherm_Model_Without_Neighbor_Interactions
                    Run_Adsorption_Simulation(KMCModel_Setup,
                                              TrackingArrays.temp_tau,
                                              TrackingArrays.temp_theta,
                                              TrackingArrays.temp_charge,
                                              TrackingArrays.plot_sites)

                Case Simulation_Type.Isotherm_Model_With_Neighbor_Interaction_Blocking
                    Run_Adsorption_Interaction_Simulation(KMCModel_Setup,
                                                          TrackingArrays.temp_tau,
                                                          TrackingArrays.temp_theta,
                                                          TrackingArrays.temp_charge,
                                                          TrackingArrays.plot_sites)

                Case Simulation_Type.Isotherm_Model_With_Diffusion
            End Select

            Dim charge_idx As Integer = 0
            For idx2 As Integer = 0 To TrackingArrays.temp_tau.Length - 1
                TrackingArrays.tau(rep_counter, idx2) = TrackingArrays.temp_tau(idx2)
                TrackingArrays.theta(rep_counter, idx2) = TrackingArrays.temp_theta(idx2)

                If TrackingArrays.temp_charge(idx2) >= 0.000000000000001 Then
                    TrackingArrays.charge_tau(rep_counter, charge_idx) = TrackingArrays.temp_theta(idx2)
                    TrackingArrays.charge(rep_counter, charge_idx) = (TrackingArrays.temp_charge(idx2) * charge_electron) / KMCModel_Setup.computational_box_dimensions.base_area
                    charge_idx = charge_idx + 1
                End If

            Next
            Console.Write("#")
        Next
        Console.WriteLine()
        Console.WriteLine("Now, writing output files...")
        ' Output the final configuration of the adsorbed species
        Console.WriteLine("Final adsorption configuration")
        WriteConfigurationFile(KMCModel_Setup.final_configuration_filename, TrackingArrays.plot_sites)

        'Output the desorption current density
        Console.WriteLine("Current density vs time...")
        WriteCurrentDensityFile(KMCModel_Setup.charge_transfer_filename,
                                KMCModel_Setup.total_time_steps_per_simulation,
                                KMCModel_Setup.number_of_simulation_repetitions,
                                TrackingArrays.tau,
                                TrackingArrays.charge)

        If KMCModel_Setup.analytical_model_available = True Then
            ' Compute the analytical solution
            Dim delta_analytical_tau As Double
            delta_analytical_tau = TrackingArrays.tau(KMCModel_Setup.number_of_simulation_repetitions - 1, KMCModel_Setup.total_time_steps_per_simulation - 1) /
                        KMCModel_Setup.total_time_steps_per_simulation
            TrackingArrays.analytical_theta(0) = 0.0
            TrackingArrays.analytical_tau(0) = 0.0
            For idx As Integer = 0 To KMCModel_Setup.total_time_steps_per_simulation - 1
                TrackingArrays.analytical_tau(idx + 1) = TrackingArrays.analytical_tau(idx) + delta_analytical_tau
                TrackingArrays.analytical_theta(idx + 1) = (KMCModel_Setup.reaction_list(0).k_reaction / (KMCModel_Setup.reaction_list(0).k_reaction + KMCModel_Setup.reaction_list(1).k_reaction)) *
                            (1 - Math.Exp(-(KMCModel_Setup.reaction_list(0).k_reaction + KMCModel_Setup.reaction_list(1).k_reaction) * TrackingArrays.analytical_tau(idx + 1)))
            Next
            Console.WriteLine("Theta vs time with analytic solution...")
            WriteCoverageFileWithAnalyticSolution(KMCModel_Setup.coverage_filename,
                                                          KMCModel_Setup.total_time_steps_per_simulation,
                                                          KMCModel_Setup.number_of_simulation_repetitions,
                                                          TrackingArrays.analytical_tau,
                                                          TrackingArrays.analytical_theta,
                                                          TrackingArrays.tau,
                                                          TrackingArrays.theta)
        Else
            Console.WriteLine("Theta vs time without analytic solution...")
            WriteCoverageFileWithoutAnalyticSolution(KMCModel_Setup.coverage_filename,
                                                             KMCModel_Setup.total_time_steps_per_simulation,
                                                             KMCModel_Setup.number_of_simulation_repetitions,
                                                             TrackingArrays.tau,
                                                             TrackingArrays.theta)
        End If

        Console.WriteLine("Finished!")
        Console.WriteLine("Press any key to exit.")
        Console.ReadKey()


    End Sub

    Sub BuildReactionLattice(ByVal simulation As ModelProperties,
                        ByRef lattice() As ReactionPosition)

        Dim lattice_site_counter As Integer = 0
        Dim total_lattice_sites, nx, ny As Integer

        'This set of nested loops builds the lattice
        ' ===========================================
        For z_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_z(0) - 1
            For y_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_y(0) - 1
                For x_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_x(0) - 1

                    If simulation.computational_box_dimensions.num_cells_z(0) > 1 Then
                        If z_counter < simulation.computational_box_dimensions.num_cells_z(0) - 1 Then
                            For b_counter As Integer = 0 To simulation.reaction_surface.basis - 1
                                lattice(lattice_site_counter) = New ReactionPosition With {
                                    .Index = lattice_site_counter
                                }

                                lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                                (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0))

                                lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                                (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) 'lattice_parameter(1)

                                lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                                (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) 'lattice_parameter(2)

                                lattice(lattice_site_counter).occupation_state = Occupation_Type.Metal_Site_Metal_Atom

                                lattice_site_counter += 1

                                For ad_site_counter As Integer = 0 To simulation.reaction_surface.num_adsorption_sites - 1
                                    lattice(lattice_site_counter) = New ReactionPosition With {
                                        .Index = lattice_site_counter
                                    }

                                    lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                                        (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0)) +
                                    (simulation.reaction_surface.adsorption_basis_x(ad_site_counter) * simulation.reaction_surface.lattice_parameter(0))

                                    lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                                        (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) +
                                    (simulation.reaction_surface.adsorption_basis_y(ad_site_counter) * simulation.reaction_surface.lattice_parameter(1))

                                    lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                                    (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) +
                                    (simulation.reaction_surface.adsorption_basis_z(ad_site_counter) * simulation.reaction_surface.lattice_parameter(2))

                                    lattice(lattice_site_counter).occupation_state = Occupation_Type.Adsorption_Site_Blocked
                                    ReDim lattice(lattice_site_counter).neighbors(simulation.reaction_surface.num_adsorption_neighbors)

                                    lattice_site_counter += 1
                                Next

                            Next
                        Else
                            For b_counter As Integer = 0 To simulation.reaction_surface.basis - 1
                                lattice(lattice_site_counter) = New ReactionPosition With {
                                    .Index = lattice_site_counter
                                }

                                lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                                (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0))

                                lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                                (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) 'lattice_parameter(1)

                                lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                                (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) 'lattice_parameter(2)

                                lattice(lattice_site_counter).occupation_state = Occupation_Type.Metal_Site_Metal_Atom

                                lattice_site_counter += 1

                                For ad_site_counter As Integer = 0 To simulation.reaction_surface.num_adsorption_sites - 1
                                    lattice(lattice_site_counter) = New ReactionPosition With {
                                        .Index = lattice_site_counter
                                    }

                                    lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                                        (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0)) +
                                    (simulation.reaction_surface.adsorption_basis_x(ad_site_counter) * simulation.reaction_surface.lattice_parameter(0))

                                    lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                                        (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) +
                                    (simulation.reaction_surface.adsorption_basis_y(ad_site_counter) * simulation.reaction_surface.lattice_parameter(1))

                                    lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                                    (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) +
                                    (simulation.reaction_surface.adsorption_basis_z(ad_site_counter) * simulation.reaction_surface.lattice_parameter(2))

                                    lattice(lattice_site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty
                                    ReDim lattice(lattice_site_counter).neighbors(simulation.reaction_surface.num_adsorption_neighbors)

                                    lattice_site_counter += 1
                                Next

                            Next
                        End If

                    Else
                        For b_counter As Integer = 0 To simulation.reaction_surface.basis - 1
                            lattice(lattice_site_counter) = New ReactionPosition With {
                                .Index = lattice_site_counter
                            }

                            lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                            (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0))

                            lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                            (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) 'lattice_parameter(1)

                            lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                            (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) 'lattice_parameter(2)

                            lattice(lattice_site_counter).occupation_state = Occupation_Type.Metal_Site_Metal_Atom

                            lattice_site_counter += 1

                            For ad_site_counter As Integer = 0 To simulation.reaction_surface.num_adsorption_sites - 1
                                lattice(lattice_site_counter) = New ReactionPosition With {
                                    .Index = lattice_site_counter
                                }

                                lattice(lattice_site_counter).position(0) = (x_counter * simulation.reaction_surface.lattice_parameter(0)) +
                                    (simulation.reaction_surface.basis_vector_x(b_counter) * simulation.reaction_surface.lattice_parameter(0)) +
                                (simulation.reaction_surface.adsorption_basis_x(ad_site_counter) * simulation.reaction_surface.lattice_parameter(0))

                                lattice(lattice_site_counter).position(1) = (y_counter * simulation.reaction_surface.lattice_parameter(1)) +
                                    (simulation.reaction_surface.basis_vector_y(b_counter) * simulation.reaction_surface.lattice_parameter(1)) +
                                (simulation.reaction_surface.adsorption_basis_y(ad_site_counter) * simulation.reaction_surface.lattice_parameter(1))

                                lattice(lattice_site_counter).position(2) = (z_counter * simulation.reaction_surface.lattice_parameter(2)) +
                                (simulation.reaction_surface.basis_vector_z(b_counter) * simulation.reaction_surface.lattice_parameter(2)) +
                                (simulation.reaction_surface.adsorption_basis_z(ad_site_counter) * simulation.reaction_surface.lattice_parameter(2))

                                lattice(lattice_site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty
                                ReDim lattice(lattice_site_counter).neighbors(simulation.reaction_surface.num_adsorption_neighbors)

                                lattice_site_counter += 1
                            Next

                        Next

                    End If

                Next
            Next
        Next

        If simulation.type_of_simulation <> Simulation_Type.Isotherm_Model_Without_Neighbor_Interactions Then
            ' ===========================================
            'This set of nested loops finds the neighbors for the adsorption sites
            ' ===========================================

            total_lattice_sites = lattice_site_counter
            lattice_site_counter = 0
            nx = CInt(simulation.computational_box_dimensions.num_cells_x(0)) '/ 2
            ny = CInt(simulation.computational_box_dimensions.num_cells_y(0)) ' / 2
            For z_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_z(0) - 1
                For y_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_y(0) - 1
                    For x_counter As Integer = 0 To simulation.computational_box_dimensions.num_cells_x(0) - 1
                        For b_counter As Integer = 0 To simulation.reaction_surface.basis - 1
                            lattice_site_counter += 1
                            For ad_site_counter As Integer = 0 To simulation.reaction_surface.num_adsorption_sites - 1

                                'If lattice_site_counter = 20103 Then
                                '    Console.WriteLine(lattice(lattice_site_counter).occupation_state)
                                'End If
                                If lattice(lattice_site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty Or
                                        lattice(lattice_site_counter).occupation_state = Occupation_Type.Adsorption_Site_Blocked Then

                                    If y_counter = 0 Then
                                        If x_counter = 0 Then ' Bottom-left corner
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2) + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx) + (ny * 2 * nx)
                                        ElseIf x_counter = simulation.computational_box_dimensions.num_cells_x(0) - 1 Then ' Bottom-right corner
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2) - (2 * nx)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx) + (ny * 2 * nx)
                                        ElseIf x_counter > 0 And x_counter < simulation.computational_box_dimensions.num_cells_x(0) - 1 Then ' Bottom side
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx) + (ny * 2 * nx)
                                        End If
                                    ElseIf y_counter = simulation.computational_box_dimensions.num_cells_y(0) - 1 Then
                                        If x_counter = 0 Then 'Top-left corner
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx) - (ny * 2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2) + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        ElseIf x_counter = simulation.computational_box_dimensions.num_cells_x(0) - 1 Then 'Top-right corner
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2) - (2 * nx)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx) - (ny * 2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        ElseIf x_counter > 0 And x_counter < simulation.computational_box_dimensions.num_cells_x(0) - 1 Then 'Top side
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx) - (ny * 2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        End If
                                    ElseIf y_counter > 0 And y_counter < simulation.computational_box_dimensions.num_cells_y(0) Then
                                        If x_counter = 0 Then 'Left side
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2) + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        ElseIf x_counter = simulation.computational_box_dimensions.num_cells_x(0) - 1 Then 'Right side
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2) - (2 * nx)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        ElseIf x_counter > 0 And x_counter < simulation.computational_box_dimensions.num_cells_x(0) - 1 Then 'Main area
                                            lattice(lattice_site_counter).neighbors(0) = lattice_site_counter + (2)
                                            lattice(lattice_site_counter).neighbors(1) = lattice_site_counter + (2 * nx)
                                            lattice(lattice_site_counter).neighbors(2) = lattice_site_counter - (2)
                                            lattice(lattice_site_counter).neighbors(3) = lattice_site_counter - (2 * nx)
                                        End If

                                    End If

                                End If
                                lattice_site_counter = lattice_site_counter + 1
                            Next
                        Next
                    Next
                Next
            Next
        End If
        ' ===========================================
    End Sub

    Sub Run_Adsorption_Simulation(ByVal simulation_setup As ModelProperties,
                                      ByRef time() As Double,
                                      ByRef kmc_theta() As Double,
                                  ByRef kmc_charge() As Double,
                                      ByRef reaction_sites_to_plot As LinkedList(Of ReactionPosition))

        Dim sites_lists As Set_of_LinkedLists = New Set_of_LinkedLists With {
            .adsorption_sites_for_reactions = New LinkedList(Of ReactionPosition),
        .inert_reaction_sites = New LinkedList(Of ReactionPosition),
        .reaction_sites_blocked = New LinkedList(Of ReactionPosition),
        .reaction_sites_empty = New LinkedList(Of ReactionPosition),
        .reaction_sites_filled = New LinkedList(Of ReactionPosition),
        .only_sites_for_reactions = New LinkedList(Of ReactionPosition)}

        Dim total_reaction_sites_available, total_sites_in_simulation As Integer

        ReDim time(simulation_setup.total_time_steps_per_simulation)
        ReDim kmc_theta(simulation_setup.total_time_steps_per_simulation)
        ReDim kmc_charge(simulation_setup.total_time_steps_per_simulation)

        'Initialize the lattice
        total_sites_in_simulation = (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.basis) +
            (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.num_adsorption_sites)

        Dim lattice_position(total_sites_in_simulation) As ReactionPosition

        BuildReactionLattice(simulation_setup, lattice_position)

        'Build a list of the empty and occupied sites
        For site_counter As Integer = 0 To total_sites_in_simulation - 1
            If lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
                sites_lists.reaction_sites_empty.AddLast(lattice_position(site_counter))
                sites_lists.only_sites_for_reactions.AddLast(lattice_position(site_counter))
            ElseIf lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                sites_lists.reaction_sites_filled.AddLast(lattice_position(site_counter))
            Else
                sites_lists.inert_reaction_sites.AddLast(lattice_position(site_counter))
            End If
        Next

        total_reaction_sites_available = sites_lists.reaction_sites_empty.Count + sites_lists.reaction_sites_filled.Count
        Determine_Transition_Probabilities(simulation_setup.reaction_list,
                                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                                           simulation_setup.W_probs) 'This initializes the number of transition probabilities to be considered

        Dim total_rate, random_interval As Double

        Randomize()
        Dim test_site_index, total_sites_to_test, track_index As Integer
        Dim random_test As Double

        time(0) = 0.0
        kmc_theta(0) = 0.0
        kmc_charge(0) = 0.0

        For iteration_counter As Integer = 0 To simulation_setup.total_time_steps_per_simulation - 1
            'Compute rate weighting for time calculation
            total_rate = (sites_lists.reaction_sites_filled.Count * simulation_setup.reaction_list(1).k_reaction) + (sites_lists.reaction_sites_empty.Count * simulation_setup.reaction_list(0).k_reaction)

            ' Select a random site
            total_sites_to_test = sites_lists.reaction_sites_empty.Count + sites_lists.reaction_sites_filled.Count
            test_site_index = CInt(Int(((total_sites_to_test) * Rnd())))
            track_index = sites_lists.only_sites_for_reactions(test_site_index).index

            'Check occupancy of the test site
            If sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
                random_test = Rnd()
                If random_test <= simulation_setup.W_probs(0) Then
                    lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
                    sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
                    sites_lists.reaction_sites_filled.AddLast(lattice_position(track_index))
                    sites_lists.reaction_sites_empty.Remove(lattice_position(track_index))
                    'kmc_charge(iteration_counter + 1) = kmc_charge(iteration_counter)
                End If
            ElseIf sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                random_test = Rnd()
                If random_test <= simulation_setup.W_probs(1) Then
                    lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
                    sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
                    sites_lists.reaction_sites_empty.AddLast(lattice_position(track_index))
                    sites_lists.reaction_sites_filled.Remove(lattice_position(track_index))
                    kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).k_reaction * sites_lists.reaction_sites_filled.Count 'simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
                End If

            End If

            random_interval = Rnd()
            time(iteration_counter + 1) = time(iteration_counter) + (-1 / total_rate) * Math.Log(random_interval)
            kmc_theta(iteration_counter + 1) = sites_lists.reaction_sites_filled.Count / total_reaction_sites_available

            Determine_Transition_Probabilities(simulation_setup.reaction_list,
                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                           simulation_setup.W_probs)
        Next

        For track_counter As Integer = 0 To sites_lists.reaction_sites_filled.Count - 1
            reaction_sites_to_plot.AddLast(sites_lists.reaction_sites_filled(track_counter))
        Next
        For track_counter2 As Integer = 0 To sites_lists.inert_reaction_sites.Count - 1
            reaction_sites_to_plot.AddLast(sites_lists.inert_reaction_sites(track_counter2))
        Next


    End Sub

    Sub Run_Adsorption_Interaction_Simulation(ByVal simulation_setup As ModelProperties,
                                      ByRef time() As Double,
                                      ByRef kmc_theta() As Double,
                                        ByRef kmc_charge() As Double,
                                      ByRef reaction_sites_to_plot As LinkedList(Of ReactionPosition))

        Dim sites_lists As Set_of_LinkedLists = New Set_of_LinkedLists With {
            .adsorption_sites_for_reactions = New LinkedList(Of ReactionPosition),
        .inert_reaction_sites = New LinkedList(Of ReactionPosition),
        .reaction_sites_blocked = New LinkedList(Of ReactionPosition),
        .reaction_sites_empty = New LinkedList(Of ReactionPosition),
        .reaction_sites_filled = New LinkedList(Of ReactionPosition)}

        Dim total_reaction_sites_available, total_sites_in_simulation As Integer
        ReDim time(simulation_setup.total_time_steps_per_simulation)
        ReDim kmc_theta(simulation_setup.total_time_steps_per_simulation)

        'Initialize the lattice
        total_sites_in_simulation = (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.basis) +
            (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.num_adsorption_sites)

        Dim lattice_position(total_sites_in_simulation) As ReactionPosition

        BuildReactionLattice(simulation_setup, lattice_position)

        'Build a list of the empty and occupied sites before starting
        For site_counter As Integer = 0 To total_sites_in_simulation - 1
            If lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
                sites_lists.reaction_sites_empty.AddLast(lattice_position(site_counter))
                sites_lists.adsorption_sites_for_reactions.AddLast(lattice_position(site_counter))
            ElseIf lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                sites_lists.reaction_sites_filled.AddLast(lattice_position(site_counter))
                sites_lists.adsorption_sites_for_reactions.AddLast(lattice_position(site_counter))
            Else
                sites_lists.inert_reaction_sites.AddLast(lattice_position(site_counter))
            End If
        Next

        total_reaction_sites_available = sites_lists.adsorption_sites_for_reactions.Count 'reaction_sites_empty.Count + reaction_sites_filled.Count
        Determine_Transition_Probabilities(simulation_setup.reaction_list,
                                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                                           simulation_setup.W_probs) 'This initializes the number of transition probabilities to be considered

        Dim total_rate, random_interval As Double

        Randomize()
        Dim test_site_index, total_sites_to_test, track_index, count_of_blocked_sites As Integer
        Dim random_test As Double
        Dim neighbors_of_selected_site(simulation_setup.reaction_surface.num_adsorption_neighbors) As Integer

        time(0) = 0.0
        kmc_theta(0) = 0.0
        kmc_charge(0) = 0.0

        For iteration_counter As Integer = 0 To simulation_setup.total_time_steps_per_simulation - 1
            'Compute rate weighting for time calculation
            total_rate = (sites_lists.reaction_sites_filled.Count * simulation_setup.reaction_list(1).k_reaction) + (sites_lists.reaction_sites_empty.Count * simulation_setup.reaction_list(0).k_reaction)
            count_of_blocked_sites = sites_lists.reaction_sites_blocked.Count

            ' Select a random site
            total_sites_to_test = sites_lists.reaction_sites_empty.Count + sites_lists.reaction_sites_filled.Count ' 5 '

            test_site_index = CInt(Int(((total_sites_to_test) * Rnd())))

            track_index = sites_lists.adsorption_sites_for_reactions(test_site_index).index 'test_sites_for_reactions(test_site_index).index '
            For idx As Integer = 0 To simulation_setup.reaction_surface.num_adsorption_neighbors - 1
                neighbors_of_selected_site(idx) = lattice_position(track_index).neighbors(idx)
            Next

            'Check occupancy of the test site
            If lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
                random_test = Rnd()
                If random_test <= simulation_setup.W_probs(0) Then
                    Dim all_clear As Boolean = True

                    For neighbor_check_idx As Integer = 0 To simulation_setup.reaction_surface.num_adsorption_neighbors - 1
                        If lattice_position(neighbors_of_selected_site(neighbor_check_idx)).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                            all_clear = False
                        End If
                    Next

                    If all_clear Then
                        lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2

                        'adsorption_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
                        sites_lists.reaction_sites_filled.AddLast(lattice_position(track_index))
                        sites_lists.reaction_sites_empty.Remove(lattice_position(track_index))
                        'kmc_charge(iteration_counter + 1) = kmc_charge(iteration_counter)
                    End If

                End If

            ElseIf lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                random_test = Rnd()
                If random_test <= simulation_setup.W_probs(1) Then
                    lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty

                    'adsorption_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
                    sites_lists.reaction_sites_empty.AddLast(lattice_position(track_index))
                    sites_lists.reaction_sites_filled.Remove(lattice_position(track_index))
                    'kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
                    kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).k_reaction * sites_lists.reaction_sites_filled.Count 'simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
                End If

            End If

            For didx As Integer = 1 To sites_lists.reaction_sites_empty.Count
                sites_lists.adsorption_sites_for_reactions.AddLast(sites_lists.reaction_sites_empty(didx - 1))
            Next
            For cidx As Integer = 1 To sites_lists.reaction_sites_filled.Count
                sites_lists.adsorption_sites_for_reactions.AddLast(sites_lists.reaction_sites_filled(cidx - 1))
            Next

            random_interval = Rnd()
            time(iteration_counter + 1) = time(iteration_counter) + (-1 / total_rate) * Math.Log(random_interval)
            kmc_theta(iteration_counter + 1) = sites_lists.reaction_sites_filled.Count / total_reaction_sites_available

            Determine_Transition_Probabilities(simulation_setup.reaction_list,
                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                           simulation_setup.W_probs)
        Next


        For track_counter As Integer = 0 To sites_lists.reaction_sites_filled.Count - 1
            reaction_sites_to_plot.AddLast(sites_lists.reaction_sites_filled(track_counter))
        Next
        For track_counter2 As Integer = 0 To sites_lists.inert_reaction_sites.Count - 1
            reaction_sites_to_plot.AddLast(sites_lists.inert_reaction_sites(track_counter2))
        Next

    End Sub

    Function Determine_Reaction_Rates(ByVal temperature As Double, ByVal reaction_type As Kinetics_Types) As Reaction_Process

        Dim reaction_rate As Reaction_Process = New Reaction_Process
        Dim q_fraction As Double

        Select Case reaction_type
            Case Kinetics_Types.Adsorption_Basic
                reaction_rate.k_reaction = 1.0 ' (sites sec)^(-1)
                reaction_rate.z_electrons = 0

            Case Kinetics_Types.Adsorption_WithCharge
                reaction_rate.prefactor = (Boltzmann_Constant * temperature) / Planck_Constant
                reaction_rate.q_AB = 10.0
                reaction_rate.q_A = 20000000000.0
                reaction_rate.q_B = 20000000000.0
                q_fraction = reaction_rate.q_AB / (reaction_rate.q_A * reaction_rate.q_B)
                reaction_rate.RT = Ideal_Gas_Constant * temperature
                reaction_rate.del_H = 0.1 * reaction_rate.RT
                reaction_rate.k_reaction = reaction_rate.prefactor * q_fraction * Math.Exp(-reaction_rate.del_H / reaction_rate.RT)
                reaction_rate.z_electrons = 0

            Case Kinetics_Types.Desorption_Basic
                reaction_rate.k_reaction = 2.0 ' (sites sec)^(-1)
                reaction_rate.z_electrons = 0

            Case Kinetics_Types.Desorption_WithCharge
                reaction_rate.prefactor = (Boltzmann_Constant * temperature) / Planck_Constant
                reaction_rate.q_AB = 20.0
                reaction_rate.q_A = 20000000000.0
                reaction_rate.q_B = 20000000000.0
                q_fraction = reaction_rate.q_AB / (reaction_rate.q_A * reaction_rate.q_B)
                reaction_rate.RT = Ideal_Gas_Constant * temperature
                reaction_rate.del_H = 0.1 * reaction_rate.RT
                reaction_rate.k_reaction = reaction_rate.prefactor * q_fraction * Math.Exp(-reaction_rate.del_H / reaction_rate.RT)
                reaction_rate.z_electrons = 4

            Case Kinetics_Types.Surface_Hopping
        End Select

        Return reaction_rate
    End Function

    Sub Determine_Transition_Probabilities(ByVal reactions() As Reaction_Process, ByVal sites_free As Double, byval_sites_occupied As Double, ByRef prob_density() As Double)
        'Generate the transition probabilities
        Dim max_rate As Double

        max_rate = 0.0
        For idx As Integer = 0 To reactions.Length - 1
            If Math.Abs(reactions(idx).k_reaction) >= max_rate Then
                max_rate = reactions(idx).k_reaction
            End If
        Next

        For idx2 As Integer = 0 To reactions.Length - 1
            prob_density(idx2) = reactions(idx2).k_reaction / max_rate
        Next

    End Sub


End Module


Steven Policastro
5:15 PM (1 minute ago)
to me 

    Sub BuildElectrolyteReservoir(ByVal box_dimensions As ComputationalCellParameters,
                                  ByVal initial_concentration As Double,
                                  ByRef reduction_species_capacity_per_time As Integer)

    reduction_species_capacity_per_time = 10 'per second

End Sub

Sub Run_Adsorption_Simulation(ByVal simulation_setup As ModelProperties,
                                  ByRef time() As Double,
                                  ByRef kmc_theta() As Double,
                              ByRef kmc_charge() As Double,
                                  ByRef reaction_sites_to_plot As LinkedList(Of ReactionPosition))

    Dim sites_lists As Set_of_LinkedLists = New Set_of_LinkedLists With {
        .adsorption_sites_for_reactions = New LinkedList(Of ReactionPosition),
    .inert_reaction_sites = New LinkedList(Of ReactionPosition),
    .reaction_sites_blocked = New LinkedList(Of ReactionPosition),
    .reaction_sites_empty = New LinkedList(Of ReactionPosition),
    .reaction_sites_filled = New LinkedList(Of ReactionPosition),
    .only_sites_for_reactions = New LinkedList(Of ReactionPosition)}

    Dim total_reaction_sites_available, total_sites_in_simulation As Integer

    ReDim time(simulation_setup.total_time_steps_per_simulation)
    ReDim kmc_theta(simulation_setup.total_time_steps_per_simulation)
    ReDim kmc_charge(simulation_setup.total_time_steps_per_simulation)

    'Initialize the lattice
    total_sites_in_simulation = (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.basis) +
        (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.num_adsorption_sites)

    Dim lattice_position(total_sites_in_simulation) As ReactionPosition

    BuildReactionLattice(simulation_setup, lattice_position)
    Dim rate_o2_available, rate_o2_demanded As Integer
    Dim rate_reduction As Double = 1.0
    BuildElectrolyteReservoir(simulation_setup.computational_box_dimensions, 0.6, rate_o2_available)
    'Adsorption process
    If random_test <= simulation_setup.W_probs(0) Then
        rate_o2_demanded = CInt(simulation_setup.reaction_list(0).k_reaction * sites_lists.reaction_sites_empty.Count)

        If rate_o2_demanded <= (rate_o2_available * sites_lists.reaction_sites_empty.Count) Then
            lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
            sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
            sites_lists.reaction_sites_filled.AddLast(lattice_position(track_index))
            sites_lists.reaction_sites_empty.Remove(lattice_position(track_index))
            rate_reduction = rate_reduction - 0.001
            rate_o2_available = CInt(rate_reduction * rate_o2_available)
            'kmc_charge(iteration_counter + 1) = kmc_charge(iteration_counter)
        End If

    End If
    ElseIf sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
    random_test = Rnd()
    'Desorption process
    If random_test <= simulation_setup.W_probs(1) Then
        lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
        sites_lists.only_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
        sites_lists.reaction_sites_empty.AddLast(lattice_position(track_index))
        sites_lists.reaction_sites_filled.Remove(lattice_position(track_index))
        kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).k_reaction * sites_lists.reaction_sites_filled.Count 'simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
    End If

    End If

    random_interval = Rnd()
    time(iteration_counter + 1) = time(iteration_counter) + (-1 / total_rate) * Math.Log(random_interval)
    kmc_theta(iteration_counter + 1) = sites_lists.reaction_sites_filled.Count / total_reaction_sites_available

    Determine_Transition_Probabilities(simulation_setup.reaction_list,
                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                           simulation_setup.W_probs)
    Next

    For track_counter As Integer = 0 To sites_lists.reaction_sites_filled.Count - 1
        reaction_sites_to_plot.AddLast(sites_lists.reaction_sites_filled(track_counter))
    Next
    For track_counter2 As Integer = 0 To sites_lists.inert_reaction_sites.Count - 1
        reaction_sites_to_plot.AddLast(sites_lists.inert_reaction_sites(track_counter2))
    Next


End Sub

Sub Run_Adsorption_Interaction_Simulation(ByVal simulation_setup As ModelProperties,
                                      ByRef time() As Double,
                                      ByRef kmc_theta() As Double,
                                        ByRef kmc_charge() As Double,
                                      ByRef reaction_sites_to_plot As LinkedList(Of ReactionPosition))

    Dim sites_lists As Set_of_LinkedLists = New Set_of_LinkedLists With {
            .adsorption_sites_for_reactions = New LinkedList(Of ReactionPosition),
        .inert_reaction_sites = New LinkedList(Of ReactionPosition),
        .reaction_sites_blocked = New LinkedList(Of ReactionPosition),
        .reaction_sites_empty = New LinkedList(Of ReactionPosition),
        .reaction_sites_filled = New LinkedList(Of ReactionPosition)}

    Dim total_reaction_sites_available, total_sites_in_simulation As Integer
    ReDim time(simulation_setup.total_time_steps_per_simulation)
    ReDim kmc_theta(simulation_setup.total_time_steps_per_simulation)

    'Initialize the lattice
    total_sites_in_simulation = (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.basis) +
            (simulation_setup.computational_box_dimensions.num_cells_x(0) * simulation_setup.computational_box_dimensions.num_cells_y(0) * simulation_setup.computational_box_dimensions.num_cells_z(0) * simulation_setup.reaction_surface.num_adsorption_sites)

    Dim lattice_position(total_sites_in_simulation) As ReactionPosition

    BuildReactionLattice(simulation_setup, lattice_position)

    'Build a list of the empty and occupied sites before starting
    For site_counter As Integer = 0 To total_sites_in_simulation - 1
        If lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
            sites_lists.reaction_sites_empty.AddLast(lattice_position(site_counter))
            sites_lists.adsorption_sites_for_reactions.AddLast(lattice_position(site_counter))
        ElseIf lattice_position(site_counter).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
            sites_lists.reaction_sites_filled.AddLast(lattice_position(site_counter))
            sites_lists.adsorption_sites_for_reactions.AddLast(lattice_position(site_counter))
        Else
            sites_lists.inert_reaction_sites.AddLast(lattice_position(site_counter))
        End If
    Next

    total_reaction_sites_available = sites_lists.adsorption_sites_for_reactions.Count 'reaction_sites_empty.Count + reaction_sites_filled.Count
    Determine_Transition_Probabilities(simulation_setup.reaction_list,
                                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                                           simulation_setup.W_probs) 'This initializes the number of transition probabilities to be considered

    Dim total_rate, random_interval As Double

    Randomize()
    Dim test_site_index, total_sites_to_test, track_index, count_of_blocked_sites As Integer
    Dim random_test As Double
    Dim neighbors_of_selected_site(simulation_setup.reaction_surface.num_adsorption_neighbors) As Integer

    time(0) = 0.0
    kmc_theta(0) = 0.0
    kmc_charge(0) = 0.0

    For iteration_counter As Integer = 0 To simulation_setup.total_time_steps_per_simulation - 1
        'Compute rate weighting for time calculation
        total_rate = (sites_lists.reaction_sites_filled.Count * simulation_setup.reaction_list(1).k_reaction) + (sites_lists.reaction_sites_empty.Count * simulation_setup.reaction_list(0).k_reaction)
        count_of_blocked_sites = sites_lists.reaction_sites_blocked.Count

        ' Select a random site
        total_sites_to_test = sites_lists.reaction_sites_empty.Count + sites_lists.reaction_sites_filled.Count ' 5 '

        test_site_index = CInt(Int(((total_sites_to_test) * Rnd())))

        track_index = sites_lists.adsorption_sites_for_reactions(test_site_index).index 'test_sites_for_reactions(test_site_index).index '
        For idx As Integer = 0 To simulation_setup.reaction_surface.num_adsorption_neighbors - 1
            neighbors_of_selected_site(idx) = lattice_position(track_index).neighbors(idx)
        Next

        'Check occupancy of the test site
        If lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty Then
            random_test = Rnd()
            If random_test <= simulation_setup.W_probs(0) Then
                Dim all_clear As Boolean = True

                For neighbor_check_idx As Integer = 0 To simulation_setup.reaction_surface.num_adsorption_neighbors - 1
                    If lattice_position(neighbors_of_selected_site(neighbor_check_idx)).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                        all_clear = False
                    End If
                Next

                If all_clear Then
                    lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2

                    'adsorption_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2
                    sites_lists.reaction_sites_filled.AddLast(lattice_position(track_index))
                    sites_lists.reaction_sites_empty.Remove(lattice_position(track_index))
                    'kmc_charge(iteration_counter + 1) = kmc_charge(iteration_counter)
                End If

            End If

        ElseIf lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
            random_test = Rnd()
            If random_test <= simulation_setup.W_probs(1) Then
                lattice_position(track_index).occupation_state = Occupation_Type.Adsorption_Site_Empty

                'adsorption_sites_for_reactions(test_site_index).occupation_state = Occupation_Type.Adsorption_Site_Empty
                sites_lists.reaction_sites_empty.AddLast(lattice_position(track_index))
                sites_lists.reaction_sites_filled.Remove(lattice_position(track_index))
                'kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
                kmc_charge(iteration_counter + 1) = simulation_setup.reaction_list(1).k_reaction * sites_lists.reaction_sites_filled.Count 'simulation_setup.reaction_list(1).z_electrons 'kmc_charge(iteration_counter) '+ ((simulation_setup.reaction_list(1).z_electrons * charge_electron) / simulation_setup.computational_box_dimensions.base_area)
            End If

        End If

        For didx As Integer = 1 To sites_lists.reaction_sites_empty.Count
            sites_lists.adsorption_sites_for_reactions.AddLast(sites_lists.reaction_sites_empty(didx - 1))
        Next
        For cidx As Integer = 1 To sites_lists.reaction_sites_filled.Count
            sites_lists.adsorption_sites_for_reactions.AddLast(sites_lists.reaction_sites_filled(cidx - 1))
        Next

        random_interval = Rnd()
        time(iteration_counter + 1) = time(iteration_counter) + (-1 / total_rate) * Math.Log(random_interval)
        kmc_theta(iteration_counter + 1) = sites_lists.reaction_sites_filled.Count / total_reaction_sites_available

        Determine_Transition_Probabilities(simulation_setup.reaction_list,
                           sites_lists.reaction_sites_empty.Count / total_reaction_sites_available,
                           sites_lists.reaction_sites_filled.Count / total_reaction_sites_available,
                           simulation_setup.W_probs)
    Next


    For track_counter As Integer = 0 To sites_lists.reaction_sites_filled.Count - 1
        reaction_sites_to_plot.AddLast(sites_lists.reaction_sites_filled(track_counter))
    Next
    For track_counter2 As Integer = 0 To sites_lists.inert_reaction_sites.Count - 1
        reaction_sites_to_plot.AddLast(sites_lists.inert_reaction_sites(track_counter2))
    Next

End Sub

Function Determine_Reaction_Rates(ByVal temperature As Double, ByVal reaction_type As Kinetics_Types) As Reaction_Process

    Dim reaction_rate As Reaction_Process = New Reaction_Process
    Dim q_fraction As Double

    Select Case reaction_type
        Case Kinetics_Types.Adsorption_Basic
            reaction_rate.k_reaction = 1.0 ' (sites sec)^(-1)
            reaction_rate.z_electrons = 0

        Case Kinetics_Types.Adsorption_WithCharge
            reaction_rate.prefactor = (Boltzmann_Constant * temperature) / Planck_Constant
            reaction_rate.q_AB = 10.0
            reaction_rate.q_A = 20000000000.0
            reaction_rate.q_B = 20000000000.0
            q_fraction = reaction_rate.q_AB / (reaction_rate.q_A * reaction_rate.q_B)
            reaction_rate.RT = Ideal_Gas_Constant * temperature
            reaction_rate.del_H = 0.1 * reaction_rate.RT
            reaction_rate.k_reaction = reaction_rate.prefactor * q_fraction * Math.Exp(-reaction_rate.del_H / reaction_rate.RT)
            reaction_rate.z_electrons = 0

        Case Kinetics_Types.Desorption_Basic
            reaction_rate.k_reaction = 2.0 ' (sites sec)^(-1)
            reaction_rate.z_electrons = 0

        Case Kinetics_Types.Desorption_WithCharge
            reaction_rate.prefactor = (Boltzmann_Constant * temperature) / Planck_Constant
            reaction_rate.q_AB = 20.0
            reaction_rate.q_A = 20000000000.0
            reaction_rate.q_B = 20000000000.0
            q_fraction = reaction_rate.q_AB / (reaction_rate.q_A * reaction_rate.q_B)
            reaction_rate.RT = Ideal_Gas_Constant * temperature
            reaction_rate.del_H = 0.1 * reaction_rate.RT
            reaction_rate.k_reaction = reaction_rate.prefactor * q_fraction * Math.Exp(-reaction_rate.del_H / reaction_rate.RT)
            reaction_rate.z_electrons = 4

        Case Kinetics_Types.Surface_Hopping
    End Select

    Return reaction_rate
End Function

Sub Determine_Transition_Probabilities(ByVal reactions() As Reaction_Process, ByVal sites_free As Double, byval_sites_occupied As Double, ByRef prob_density() As Double)
    'Generate the transition probabilities
    Dim max_rate As Double

    max_rate = 0.0
    For idx As Integer = 0 To reactions.Length - 1
        If Math.Abs(reactions(idx).k_reaction) >= max_rate Then
            max_rate = reactions(idx).k_reaction
        End If
    Next

    For idx2 As Integer = 0 To reactions.Length - 1
        prob_density(idx2) = reactions(idx2).k_reaction / max_rate
    Next

End Sub


End Module
