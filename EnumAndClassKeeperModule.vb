Imports System.IO

Module EnumAndClassKeeperModule
    Public Boltzmann_Constant As Double = 1.38064852E-23 'm2 kg s-2 K-1
    Public Planck_Constant As Double = 6.62607004E-34 'm2 kg / s
    Public Ideal_Gas_Constant As Double = 8.314 'J/mol K
    Public Faraday_Constant As Double = 96485.3329 'A/mol
    Public Avogadro_Number As Double = 6.02E+23 'Particles/mol
    Public charge_electron As Double = 1.6E-19 'C

    Public Enum Simulation_Type
        Isotherm_Model_Without_Neighbor_Interactions
        Isotherm_Model_With_Neighbor_Interaction_Blocking
        Isotherm_Model_With_Diffusion
    End Enum

    Public Enum Occupation_Type
        Adsorption_Site_Empty
        Adsorption_Site_Blocked
        Adsorption_Site_Adsorbed_O2
        Metal_Site_Metal_Atom
        Non_Reacting
        Bottom
        O_star
        OH_minus
        OH_star
    End Enum

    Public Enum Cell_Geometry
        Flat
        Cube
        Rectangular_Parallelpiped
        Island
        Trench
        Pit
    End Enum

    Public Enum Alloy
        Aluminum
        Iron
        Gold
    End Enum

    Public Enum Lattice_Type
        FCC
        BCC
        SC
    End Enum

    Public Enum Electrolyte
        Water
        NaCl_06M
        NaCl_3M
        NaCl_5M
    End Enum

    Public Enum Kinetics_Types
        Adsorption_Basic
        Desorption_Basic
        Adsorption_WithCharge
        Desorption_WithCharge
        Surface_Hopping
        SurfKin
        Corrosion
    End Enum

    Public Structure Reaction_Process
        Public prefactor As Double
        Public q_AB As Double
        Public q_A As Double
        Public q_B As Double
        Public del_H As Double
        Public RT As Double
        Public k_reaction As Double
        Public z_electrons As Integer
    End Structure

    Public Structure ComputationalCellParameters
        Public num_cells_x(), num_cells_y(), num_cells_z() As Integer
        Public metal_volume, electrolyte_volume, base_area As Double
    End Structure

    Public Class TrackingParameters
        Public temp_tau(), temp_theta(), temp_charge(), analytical_tau(), analytical_theta() As Double
        Public tau(,), theta(,), charge_tau(,), charge(,) As Double
        Public plot_sites As LinkedList(Of ReactionPosition)
    End Class

    Public Structure Set_of_LinkedLists
        Public reaction_sites_empty As LinkedList(Of ReactionPosition)
        Public reaction_sites_filled As LinkedList(Of ReactionPosition)
        Public reaction_sites_blocked As LinkedList(Of ReactionPosition)
        Public adsorption_sites_for_reactions As LinkedList(Of ReactionPosition)
        Public inert_reaction_sites As LinkedList(Of ReactionPosition)
        Public only_sites_for_reactions As LinkedList(Of ReactionPosition)
    End Structure

    Public Class ReactionPosition
        Public index As Integer
        Public position(2) As Double
        Public occupation_state As Occupation_Type
        Public neighbors() As Integer
    End Class

    Public Class Lattice_Properties
        Public lattice_name As Lattice_Type
        Public lattice_parameter(2) As Double
        Public basis As Integer
        Public basis_vector_x() As Double
        Public basis_vector_y() As Double
        Public basis_vector_z() As Double
        Public num_neighbors As Integer
        Public num_adsorption_sites As Integer
        Public adsorption_basis_x() As Double
        Public adsorption_basis_y() As Double
        Public adsorption_basis_z() As Double
        Public num_adsorption_neighbors As Integer
    End Class

End Module
