using ConstructiveSolidGeometry

inner_cyl = InfCylinder(Coord(0.0,0.0,0.0), Coord(0.0,0.0,1.0), 1.0)
outer_cyl = InfCylinder(Coord(0.0,0.0,0.0), Coord(0.0,0.0,1.0), 1.5)
left_wall = Plane(Coord(-2.0,0.0,0.0), Coord(1.0,0.0,0.0))
right_wall = Plane(Coord(2.0,0.0,0.0), Coord(-1.0,0.0,0.0))
top_wall = Plane(Coord(0.0,2.0,0.0), Coord(0.0,-1.0,0.0))
bottom_wall = Plane(Coord(0.0,-2.0,0.0), Coord(0.0,1.0,0.0))
up_wall = Plane(Coord(0.0,0.0,2.0), Coord(0.0,0.0,-1.0))
down_wall = Plane(Coord(0.0,0.0,-2.0), Coord(0.0,0.0,1.0))


inner_cell = Cell([Region(inner_cyl, -1)], :(1 ^ 1))
ring_cell = Cell([Region(inner_cyl, +1), Region(outer_cyl, -1)], :(1 ^ 2))
box_cell = Cell([Region(outer_cyl, +1), 
                 Region(left_wall, +1),
                 Region(right_wall, +1),
                 Region(top_wall, +1),
                 Region(bottom_wall, +1),
                 Region(up_wall, +1),
                 Region(down_wall, +1)], :(1 ^ 2 ^ 3 ^ 4 ^ 5 ^ 6 ^ 7))
bounding_box = Box(Coord(-2.0,-2.0,-2.0), Coord(2.0,2.0,2.0))
geom = Geometry([inner_cell, ring_cell, box_cell], bounding_box)

new_cells = [subdivide_cell(inner_cell,0,12,0); [ring_cell]; subdivide_cell(box_cell, 0,8,0)]
print(new_cells[1])
new_geom = Geometry(new_cells, bounding_box)

plot_geometry_2D(geom, Box(Coord(-2.0,-2.0,0.0), Coord(2.0,2.0,0.0)), 300)
plot_geometry_2D(new_geom, Box(Coord(-2.0,-2.0,0.0), Coord(2.0,2.0,0.0)), 300)