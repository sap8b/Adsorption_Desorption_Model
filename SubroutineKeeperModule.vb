Imports System.IO

Module SubroutineKeeperModule

    Sub WriteCoverageFileWithAnalyticSolution(ByVal filename As String,
                      ByVal iterations As Integer,
                      ByVal repetitions As Integer,
                      ByVal analytic_time() As Double,
                      ByVal analytic_theta() As Double,
                      ByVal kmc_tau(,) As Double,
                      ByVal kmc_theta(,) As Double)

        If File.Exists(filename) Then
            File.Delete(filename)
        End If

        Dim fS As FileStream
        fS = New FileStream(filename, FileMode.Create)
        Dim sW As StreamWriter = New StreamWriter(fS)

        Dim data_line, data_line_out, data_line_in As String
        data_line_in = ""
        For line_count As Integer = 0 To iterations - 1
            data_line_out = CStr(line_count) + " " +
            CStr(analytic_time(line_count)) + " " +
            CStr(analytic_theta(line_count)) + " "

            For rep_count As Integer = 0 To repetitions - 1
                data_line_in = data_line_in + CStr(kmc_tau(rep_count, line_count)) + " " + CStr(kmc_theta(rep_count, line_count)) + " "
            Next

            data_line = data_line_out + data_line_in
            sW.WriteLine(data_line)

            data_line_in = ""
        Next

        sW.Close()
        fS.Close()
    End Sub

    Sub WriteCurrentDensityFile(ByVal filename As String,
                                ByVal total_iterations As Integer,
                                ByVal repetitions As Integer,
                                ByVal kmc_tau(,) As Double,
                                ByVal kmc_charge(,) As Double)
        If File.Exists(filename) Then
            File.Delete(filename)
        End If

        Dim fS As FileStream
        fS = New FileStream(filename, FileMode.Create)
        Dim sW As StreamWriter = New StreamWriter(fS)

        Dim data_line, data_line_out, data_line_in As String
        data_line_in = ""
        For line_count As Integer = 0 To total_iterations - 1
            data_line_out = CStr(line_count) + " "

            For rep_count As Integer = 0 To repetitions - 1
                data_line_in = data_line_in + CStr(kmc_tau(rep_count, line_count)) + " " + CStr(kmc_charge(rep_count, line_count)) + " "
            Next

            data_line = data_line_out + data_line_in
            sW.WriteLine(data_line)

            data_line_in = ""
        Next
        sW.Close()
        fS.Close()
    End Sub

    Sub WriteCoverageFileWithoutAnalyticSolution(ByVal filename As String,
                      ByVal iterations As Integer,
                      ByVal repetitions As Integer,
                      ByVal kmc_tau(,) As Double,
                      ByVal kmc_theta(,) As Double)

        If File.Exists(filename) Then
            File.Delete(filename)
        End If

        Dim fS As FileStream
        fS = New FileStream(filename, FileMode.Create)
        Dim sW As StreamWriter = New StreamWriter(fS)

        Dim data_line, data_line_out, data_line_in As String
        data_line_in = ""
        For line_count As Integer = 0 To iterations - 1
            data_line_out = CStr(line_count) + " "
            For rep_count As Integer = 0 To repetitions - 1
                data_line_in = data_line_in + CStr(kmc_tau(rep_count, line_count)) + " " + CStr(kmc_theta(rep_count, line_count)) + " "
            Next

            data_line = data_line_out + data_line_in
            sW.WriteLine(data_line)

            data_line_in = ""
        Next

        sW.Close()
        fS.Close()
    End Sub

    Sub WriteConfigurationFile(ByVal filename As String,
                               ByVal sites As LinkedList(Of ReactionPosition))

        If File.Exists(filename) Then
            File.Delete(filename)
        End If

        Dim fS As FileStream
        fS = New FileStream(filename, FileMode.Create)
        Dim sW As StreamWriter = New StreamWriter(fS)

        Dim output_string, color_string As String

        For idx As Integer = 0 To sites.Count - 1
            If sites(idx).occupation_state = Occupation_Type.Adsorption_Site_Adsorbed_O2 Then
                color_string = "1"
            ElseIf sites(idx).occupation_state = Occupation_Type.Metal_Site_Metal_Atom Then
                color_string = "9"
            Else
                color_string = "0"
            End If

            output_string = CStr(sites(idx).position(0)) + " " +
                CStr(sites(idx).position(1)) + " " +
                CStr(sites(idx).position(2)) + " " +
                color_string

            sW.WriteLine(output_string)
        Next

        sW.Close()
        fS.Close()

    End Sub

    Sub WriteReactionRates(ByVal filename As String,
                      ByVal reaction_rates() As Reaction_Process)

        If File.Exists(filename) Then
            File.Delete(filename)
        End If

        Dim fS As FileStream
        fS = New FileStream(filename, FileMode.Create)
        Dim sW As StreamWriter = New StreamWriter(fS)

        Dim data_line, data_line_out, data_line_in As String
        data_line_out = ""
        For line_count As Integer = 0 To reaction_rates.Length - 1
            data_line_in = CStr(line_count)
            data_line_out = data_line_in + " " + CStr(reaction_rates(line_count).k_reaction)
            sW.WriteLine(data_line_out)
        Next

        sW.Close()
        fS.Close()
    End Sub

    Public Function Determine_BoxSize(ByVal geometry As Cell_Geometry, ByVal lattice_values As Lattice_Properties) As ComputationalCellParameters

        Dim x_length, y_length, z_length As Double

        Select Case geometry
            Case Cell_Geometry.Flat
                Dim Box_Dims As ComputationalCellParameters = New ComputationalCellParameters With {
            .num_cells_x = New Integer() {50},
            .num_cells_y = New Integer() {50},
            .num_cells_z = New Integer() {1}}

                x_length = Box_Dims.num_cells_x(0) * lattice_values.lattice_parameter(0)
                y_length = Box_Dims.num_cells_y(0) * lattice_values.lattice_parameter(1)
                z_length = Box_Dims.num_cells_z(0) * lattice_values.lattice_parameter(2)
                Box_Dims.metal_volume = x_length * y_length * z_length
                Box_Dims.electrolyte_volume = 0.0
                Box_Dims.base_area = x_length * y_length

                Return Box_Dims

            Case Cell_Geometry.Cube
                Dim Box_Dims As ComputationalCellParameters = New ComputationalCellParameters With {
            .num_cells_x = New Integer() {50},
            .num_cells_y = New Integer() {50},
            .num_cells_z = New Integer() {50}}

                x_length = Box_Dims.num_cells_x(0) * lattice_values.lattice_parameter(0)
                y_length = Box_Dims.num_cells_y(0) * lattice_values.lattice_parameter(1)
                z_length = Box_Dims.num_cells_z(0) * lattice_values.lattice_parameter(2)
                Box_Dims.metal_volume = x_length * y_length * z_length
                Box_Dims.electrolyte_volume = 0.0
                Box_Dims.base_area = x_length * y_length

                Return Box_Dims

            Case Cell_Geometry.Rectangular_Parallelpiped
                Dim Box_Dims As ComputationalCellParameters = New ComputationalCellParameters With {
            .num_cells_x = New Integer() {25},
            .num_cells_y = New Integer() {25},
            .num_cells_z = New Integer() {5}}

                x_length = Box_Dims.num_cells_x(0) * lattice_values.lattice_parameter(0)
                y_length = Box_Dims.num_cells_y(0) * lattice_values.lattice_parameter(1)
                z_length = Box_Dims.num_cells_z(0) * lattice_values.lattice_parameter(2)
                Box_Dims.metal_volume = x_length * y_length * z_length
                Box_Dims.electrolyte_volume = 0.0
                Box_Dims.base_area = x_length * y_length

                Return Box_Dims

            Case Else
                Console.WriteLine("Using non-specified cell geometry!")
                Dim Box_Dims As ComputationalCellParameters = New ComputationalCellParameters With {
            .num_cells_x = New Integer() {10},
            .num_cells_y = New Integer() {10},
            .num_cells_z = New Integer() {10}}

                x_length = Box_Dims.num_cells_x(0) * lattice_values.lattice_parameter(0)
                y_length = Box_Dims.num_cells_y(0) * lattice_values.lattice_parameter(1)
                z_length = Box_Dims.num_cells_z(0) * lattice_values.lattice_parameter(2)
                Box_Dims.metal_volume = x_length * y_length * z_length
                Box_Dims.electrolyte_volume = 0.0
                Box_Dims.base_area = x_length * y_length

                Return Box_Dims
        End Select



    End Function

    Public Function DetermineLatticeProperties(ByVal Element As Alloy) As Lattice_Properties
        Dim alloy_properties As Lattice_Properties = New Lattice_Properties

        Select Case Element
            Case Alloy.Aluminum
                alloy_properties.lattice_name = Lattice_Type.FCC

                alloy_properties.lattice_parameter(0) = 0.00000000040495 'm
                alloy_properties.lattice_parameter(1) = 0.00000000040495 'm
                alloy_properties.lattice_parameter(2) = 0.00000000040495 'm

                alloy_properties.num_neighbors = 12

                alloy_properties.basis = 4

                ReDim alloy_properties.basis_vector_x(alloy_properties.basis), alloy_properties.basis_vector_y(alloy_properties.basis),
                    alloy_properties.basis_vector_z(alloy_properties.basis)

                alloy_properties.basis_vector_x(0) = 0.0
                alloy_properties.basis_vector_y(0) = 0.0
                alloy_properties.basis_vector_z(0) = 0.0

                alloy_properties.basis_vector_x(1) = 0.5
                alloy_properties.basis_vector_y(1) = 0.5
                alloy_properties.basis_vector_z(1) = 0.0

                alloy_properties.basis_vector_x(2) = 0.5
                alloy_properties.basis_vector_y(2) = 0.0
                alloy_properties.basis_vector_z(2) = 0.5

                alloy_properties.basis_vector_x(3) = 0.0
                alloy_properties.basis_vector_y(3) = 0.5
                alloy_properties.basis_vector_z(3) = 0.5

                alloy_properties.num_adsorption_sites = 4

                ReDim alloy_properties.adsorption_basis_x(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_y(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_z(alloy_properties.num_adsorption_sites)

                alloy_properties.adsorption_basis_x(0) = 0.0
                alloy_properties.adsorption_basis_y(0) = 0.0
                alloy_properties.adsorption_basis_z(0) = 0.25

                alloy_properties.adsorption_basis_x(1) = 0.0
                alloy_properties.adsorption_basis_y(1) = 0.0
                alloy_properties.adsorption_basis_z(1) = 0.25

                alloy_properties.adsorption_basis_x(2) = 0.0
                alloy_properties.adsorption_basis_y(2) = 0.0
                alloy_properties.adsorption_basis_z(2) = 0.25

                alloy_properties.adsorption_basis_x(3) = 0.0
                alloy_properties.adsorption_basis_y(3) = 0.0
                alloy_properties.adsorption_basis_z(3) = 0.25

                alloy_properties.num_adsorption_neighbors = 4

            Case Alloy.Iron
                alloy_properties.lattice_name = Lattice_Type.BCC

                alloy_properties.lattice_parameter(0) = 0.00000000028665 'm
                alloy_properties.lattice_parameter(1) = 0.00000000028665 'm
                alloy_properties.lattice_parameter(2) = 0.00000000028665 'm

                alloy_properties.num_neighbors = 8

                alloy_properties.basis = 2

                ReDim alloy_properties.basis_vector_x(alloy_properties.basis), alloy_properties.basis_vector_y(alloy_properties.basis),
                    alloy_properties.basis_vector_z(alloy_properties.basis)

                alloy_properties.basis_vector_x(0) = 0.0
                alloy_properties.basis_vector_y(0) = 0.0
                alloy_properties.basis_vector_z(0) = 0.0

                alloy_properties.basis_vector_x(1) = 0.5
                alloy_properties.basis_vector_y(1) = 0.5
                alloy_properties.basis_vector_z(1) = 0.5

                alloy_properties.num_adsorption_sites = 2

                ReDim alloy_properties.adsorption_basis_x(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_y(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_z(alloy_properties.num_adsorption_sites)

                alloy_properties.adsorption_basis_x(0) = 0.0
                alloy_properties.adsorption_basis_y(0) = 0.0
                alloy_properties.adsorption_basis_z(0) = 0.5

                alloy_properties.adsorption_basis_x(1) = 0.0
                alloy_properties.adsorption_basis_y(1) = 0.0
                alloy_properties.adsorption_basis_z(1) = 0.5

                alloy_properties.num_adsorption_neighbors = 4

            Case Alloy.Gold
                alloy_properties.lattice_name = Lattice_Type.SC
                alloy_properties.lattice_parameter(0) = 0.00000000040782 'm
                alloy_properties.lattice_parameter(1) = 0.00000000040782 'm
                alloy_properties.lattice_parameter(2) = 0.00000000040782 'm
                alloy_properties.num_neighbors = 4

                alloy_properties.basis = 1

                ReDim alloy_properties.basis_vector_x(alloy_properties.basis), alloy_properties.basis_vector_y(alloy_properties.basis),
                    alloy_properties.basis_vector_z(alloy_properties.basis)

                alloy_properties.basis_vector_x(0) = 0.0
                alloy_properties.basis_vector_y(0) = 0.0
                alloy_properties.basis_vector_z(0) = 0.0

                alloy_properties.num_adsorption_sites = 1

                ReDim alloy_properties.adsorption_basis_x(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_y(alloy_properties.num_adsorption_sites),
                    alloy_properties.adsorption_basis_z(alloy_properties.num_adsorption_sites)

                alloy_properties.adsorption_basis_x(0) = 0.0
                alloy_properties.adsorption_basis_y(0) = 0.0
                alloy_properties.adsorption_basis_z(0) = 0.5

                alloy_properties.num_adsorption_neighbors = 4

        End Select

        Return alloy_properties

    End Function


End Module

