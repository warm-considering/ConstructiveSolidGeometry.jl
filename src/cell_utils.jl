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
    subdivide_cell(cell::Cell, rmeshes::Int64, θmeshes::Int64, zmeshes::Int64)

Creates a set of new cells by splitting a given cell into given numbers of equally spaced radial, angular, and axial subdivisions.

# Arguments
* `cell::Cell`: the geometry we want to plot
* `rmeshes::Int64`: The view box is an axis aligned box that defines where the picture will be taken, with both min and max z dimensions indicating the single z elevation the slice is taken at.
* `θmeshes::Int64`: The dimension is the number of pixels along the x and y axis to use, which determines the resolution of the picture.
* `zmeshes::Int64`: The index of the cell we wish to view
"""
function subdivide_cell(cell::Cell, rmeshes::Int64, θmeshes::Int64, zmeshes::Int64)
	@assert θmeshes % 2 == 0
	n_planes = Int64(θmeshes//2)
	angles = [i*pi/n_planes for i in 0:(n_planes-1)]

	θplanes = [Plane(Coord(0.0,0.0,0.0), Coord(cos(θ), sin(θ), 0.0)) for θ in angles]

	new_cells = Cell[]

	for i in eachindex(θplanes)
		if i == length(θplanes)
			plane1 = θplanes[i]
			plane2 = θplanes[1]
		else
			plane1 = θplanes[i]
			plane2 = θplanes[i+1]
		end

		plane1_ind = length(cell.regions) + 1
		plane2_ind = plane1_ind + 1

		#case 1: + plane1 ^ - plane2
		new_regions = push!(deepcopy(cell.regions), Region(plane1, +1), Region(plane2, -1))
		new_cell1 = Cell(new_regions, :(($(cell.definition) ^ $plane1_ind) ^ $plane2_ind))

		#case 2: - plane1 ^ +plane2
		new_regions = push!(deepcopy(cell.regions), Region(plane1, -1), Region(plane2, +1))
		new_cell2 = Cell(new_regions, :(($(cell.definition) ^ $plane1_ind) ^ $plane2_ind))

		push!(new_cells, new_cell1, new_cell2)
	end

	return new_cells
end