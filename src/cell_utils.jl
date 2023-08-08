"""
    is_in_cell(p::Coord, cell::Cell)

Determines if a point (such as a Ray origin) is inside a given cell
"""
function is_in_cell(p::Coord, cell::Cell)
    result = navigate_tree(p, cell.regions, cell.definition)
    return result
end

function navigate_tree(p::Coord, r::Array{Region}, ex::Expr)
    global _p = Coord(p.x, p.y, p.z)

	# Check if Complement
	if ex.args[1] == :~
		if typeof(ex.args[2]) == typeof(1)
			return ~ r[ex.args[2]]
		else
			return ~ navigate_tree(p, r, ex.args[2])
		end
	end

	if typeof(ex.args[2]) == typeof(1)
		# Case 1 - Both operands are leaves
		if typeof(ex.args[3]) == typeof(1)
			if ex.args[1] == :^
            	return r[ex.args[2]] ^ r[ex.args[3]]
			end
			if ex.args[1] == :|
            	return r[ex.args[2]] | r[ex.args[3]]
			end
		end
		# Case 2 - Left operand is leaf, right is not
		if typeof(ex.args[3]) != typeof(1)
			if ex.args[1] == :^
            	return r[ex.args[2]] ^ navigate_tree(p, r, ex.args[3])
			end
			if ex.args[1] == :|
            	return r[ex.args[2]] | navigate_tree(p, r, ex.args[3])
			end
		end
	end

	if typeof(ex.args[2]) != typeof(1)
		# Case 3 - left operand is not leaf, but right is
		if typeof(ex.args[3]) == typeof(1)
			if ex.args[1] == :^
            	return navigate_tree(p, r, ex.args[2]) ^ r[ex.args[3]]
			end
			if ex.args[1] == :|
            	return navigate_tree(p, r, ex.args[2]) | r[ex.args[3]]
			end
		end
		# Case 4 - Neither operand is a leaf
		if typeof(ex.args[3]) != typeof(1)
			if ex.args[1] == :^
            	return navigate_tree(p, r, ex.args[2]) ^ navigate_tree(p, r, ex.args[3])
			end
			if ex.args[1] == :|
            	return navigate_tree(p, r, ex.args[2]) | navigate_tree(p, r, ex.args[3])
			end
		end
	end
end

"""
    find_cell_id(p::Coord, geometry::Geometry)

Finds the cell id that a point resides within
"""
function find_cell_id(p::Coord, geometry::Geometry)
    for i = 1:length(geometry.cells)
        if is_in_cell(p, geometry.cells[i]) == true
            return i
        end
    end
    return -1
end

"""
    planar_subdivide(cell::Cell, origin::Coord, direction::Coord, meshes::Array{Float64})

Creates a set of new cells by splitting a given cell using parallel planes.

# Arguments
* `cell::Cell`: the cell to be subdivided
* `origin::Coord`: the position of the "zeroth" plane i.e. the first plane should be a distance 'meshes[1]' away from this position
* `direction::Coord`: the normal direction of the parallel planes
* `meshes::Array{Float64}`: list of boundaries to split the cell with

"""
function planar_subdivide(cell::Cell, origin::Coord, direction::Coord, meshes::Array{Float64})
	unit_direction = unitize(direction)
	planes = [Plane(origin + d*unit_direction, unit_direction) for d in meshes]

	new_cells = Cell[]

	for i in eachindex(planes)
		plane1 = planes[i]
		plane1_idx = length(cell.regions) + 1

		if i == 1
			#case 1: distance up to first plane
			new_regions = push!(deepcopy(cell.regions), Region(plane1, -1))
			new_cell = Cell(new_regions, :($(cell.definition) ^ $plane1_idx))


			push!(new_cells, new_cell)

		end
		
		if length(planes) > 1 && i < length(planes)
			#case 2: inner pairs of planes
			plane2 = planes[i+1]

			plane2_idx = plane1_idx + 1

			new_regions = push!(deepcopy(cell.regions), Region(plane1, +1), Region(plane2, -1))
			new_cell = Cell(new_regions, :(($(cell.definition) ^ $plane1_idx) ^ $plane2_idx))

			push!(new_cells, new_cell)
		end

		#case 3: distance after last plane
		if i == length(planes)
			new_regions = push!(deepcopy(cell.regions), Region(plane1, +1))
			new_cell = Cell(new_regions, :($(cell.definition) ^ $plane1_idx))

			push!(new_cells, new_cell)
		end	
	end

	return new_cells
end

"""
    planar_subdivide(cells::Array{Cell}, origin::Coord, direction::Coord, meshes::Array{Float64})

Creates a set of new cells by splitting each cell in a given set of cells using parallel planes.

# Arguments
* `cells::Array{Cell}`: the cells to be subdivided
* `origin::Coord`: the position of the "zeroth" plane i.e. the first plane should be a distance 'meshes[1]' away from this position
* `direction::Coord`: the normal direction of the parallel planes
* `meshes::Array{Float64}`: list of boundaries to split the cell with
"""
function planar_subdivide(cells::Array{Cell}, origin::Coord, direction::Coord, meshes::Array{Float64})
	new_cells = Cell[]

	for cell in cells
		append!(new_cells, planar_subdivide(cell, origin, direction, meshes))
	end

	return new_cells
end

"""
    azimuthal_subdivide(cell::Cell, pole::Coord, θmeshes::Array{Float64})

Creates a set of new cells by splitting a given cell into azimuthal sectors.

# Arguments
* `cell::Cell`: the cell to be subdivided
* `pole::Coord`: the vector around which to create the azimuthal subdivisions
* `θmeshes::Array{Float64}`: list of boundaries to split the cell with (values 0.0 < θ < π)
"""
function azimuthal_subdivide(cell::Cell, pole::Coord, θmeshes::Array{Float64})
	@assert maximum(θmeshes) < π

	unit_pole = unitize(pole)

	α = unit_pole.x
	β = unit_pole.y
	γ = unit_pole.z

	R = [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0] *
		[cos(β) 0.0 sin(β); 0.0 1.0 0.0; -sin(β) 0.0 cos(β)] *
		[1.0 0.0 0.0; 0.0 cos(γ) -sin(γ); 0.0 sin(γ) cos(γ)]

	#generate the normals for each splitting plane by rotating so that planes are parallel with the pole
	#TODO: extend Coords to support matrix multiplication
	transformed_normals = [R * [cos(θ); sin(θ); 0.0] for θ in θmeshes]

	θplanes = [Plane(Coord(0.0,0.0,0.0), Coord(norm[1], norm[2], norm[3])) for norm in transformed_normals]

	new_cells = Cell[]

	for i in eachindex(θplanes)
		if i == length(θplanes)
			plane1 = θplanes[i]
			plane2 = θplanes[1]
		else
			plane1 = θplanes[i]
			plane2 = θplanes[i+1]
		end

		plane1_idx = length(cell.regions) + 1
		plane2_idx = plane1_idx + 1

		#check if there is only 1 plane in the array
		if length(θplanes) == 1
			new_regions = push!(deepcopy(cell.regions), Region(plane1, +1))
			new_cell1 = Cell(new_regions, :($(cell.definition) ^ $plane1_idx))

			new_regions = push!(deepcopy(cell.regions), Region(plane1, -1))
			new_cell2 = Cell(new_regions, :($(cell.definition) ^ $plane1_idx))
		else
			#case 1: + plane1 ^ - plane2
			new_regions = push!(deepcopy(cell.regions), Region(plane1, +1), Region(plane2, -1))
			new_cell1 = Cell(new_regions, :(($(cell.definition) ^ $plane1_idx) ^ $plane2_idx))

			#case 2: - plane1 ^ +plane2
			new_regions = push!(deepcopy(cell.regions), Region(plane1, -1), Region(plane2, +1))
			new_cell2 = Cell(new_regions, :(($(cell.definition) ^ $plane1_idx) ^ $plane2_idx))
		end

		push!(new_cells, new_cell1, new_cell2)
	end

	return new_cells
end

"""
    azimuthal_subdivide(cells::Array{Cell}, pole::Coord, θmeshes::Array{Float64})

Creates a set of new cells by splitting each cell in a given set of cells into azimuthal sectors.

# Arguments
* `cells::Array{Cell}`: the cells to be subdivided
* `pole::Coord`: the vector around which to create the azimuthal subdivisions
* `θmeshes::Array{Float64}`: list of boundaries to split the cell with (values 0.0 < θ < π)
"""
function azimuthal_subdivide(cells::Array{Cell}, pole::Coord, θmeshes::Array{Float64})
	new_cells = Cell[]

	for cell in cells
		append!(new_cells, azimuthal_subdivide(cell, pole, θmeshes))
	end

	return new_cells
end

"""
    radial_subdivide(cell::Cell, pole::Coord, rmeshes::Array{Float64})

Creates a set of new cells by subdividing a given cell radially.
# Arguments
* `cell::Cell`: the cell to be subdivided
* `pole::Coord`: the direction parallel to the radial boundaries
* `rmeshes::Array{Float64}`: list of boundaries to split the cell with (values > 0.0)
"""
function radial_subdivide(cell::Cell, pole::Coord, rmeshes::Array{Float64})
	rcyls = [InfCylinder(Coord(0.0,0.0,0.0), pole, r) for r in rmeshes]

	new_cells = Cell[]

	for i in eachindex(rcyls)
		cyl1 = rcyls[i]
		cyl1_idx = length(cell.regions) + 1

		if i == 1
			#case 1: smallest cylinder inner
			new_regions = push!(deepcopy(cell.regions), Region(cyl1, -1))
			new_cell = Cell(new_regions, :($(cell.definition) ^ $cyl1_idx))


			push!(new_cells, new_cell)

		end
		
		if length(rcyls) > 1 && i < length(rcyls)
			#case 2: inner pairs of cylinders
			cyl2 = rcyls[i+1]

			cyl2_idx = cyl1_idx + 1

			new_regions = push!(deepcopy(cell.regions), Region(cyl1, +1), Region(cyl2, -1))
			new_cell = Cell(new_regions, :(($(cell.definition) ^ $cyl1_idx) ^ $cyl2_idx))

			push!(new_cells, new_cell)
		end

		#case 3: largest cylinder outer
		if i == length(rcyls)
			new_regions = push!(deepcopy(cell.regions), Region(cyl1, +1))
			new_cell = Cell(new_regions, :($(cell.definition) ^ $cyl_idx))

			push!(new_cells, new_cell)
		end
	end

	return new_cells
end

"""
    radial_subdivide(cells::Array{Cell}, pole::Coord, rmeshes::Array{Float64})

Creates a set of new cells by subdividing a given set of cells radially.
# Arguments
* `cells::Array{Cell}`: the cells to be subdivided
* `pole::Coord`: the direction parallel to the radial boundaries
* `rmeshes::Array{Float64}`: list of boundaries to split the cell with (values > 0.0)
"""
function  radial_subdivide(cells::Array{Cell}, pole::Coord, rmeshes::Array{Float64})
	new_cells = Cell[]

	for cell in cells
		append!(new_cells, radial_subdivide(cell, pole, rmeshes))
	end

	return new_cells
end
