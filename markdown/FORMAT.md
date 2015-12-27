####`.geo`

Geometry files `.geo` define metal layer elevations (z-direction), via elevations and which metal layers are connected from above and below, and the polygon descriptions for each layer. Currently, only boundary and path types of GDS2 are processed into polygons. All other types of GDS2 are ignored. 

The format is as the following


    #metal_layers
    layer_index, z_min, z_max 
    .
    .
    #vias
    layer_index, z_min, z_max, bottom_layer_index, top_layer_index
    layer_index
    {#polygons_in_a_layer
    [#vertices_in_a_polygon
    vertex_1
    vertex_2
    .
    .
    vertex_last]}
    .
    .
    

It is worth noting that metal layers and vias are indexed together, staring from 0. The square bracket part [] describes x- and y-coordinates of vertices for a polygon. There are total `#polygons_in_a_layer` such square brackets for a layer. Once all polygons in a single layer are described, we move to the next layer description enclosed by curly brackets {}. The file ends when all layers, including metal layers and vias, and all polygons in each layer are listed. The brackets are not written in output files explicitly. We put them here only for explanation purposes.

All numbers are integers with an implicit unit nm. Check if your GDS2 files are defined in different units.

Examples can be found under folder `caplet_geo/example`.



####`.caplet`

Instantiable basis function files `.caplet` describe fundamental shapes that are used to constitute instantiable basis functions. The format is as the following

    #conductors
    #shapes_for_conductor_1 #shapes_for_conductor_2 ...
    total_number_of_shapes
    [ST Inc XL XU YL YU ZL ZU Dir SDir Decay Ignore]
    .
    .

The square bracket records information for one shape. The file terminates when all shapes of `total_number_of_shapes` are listed. Each element is explained below:

`ST`: Shape Type is either 'A' for arch or 'F' for flat.

`Inc`: Increment value is either '1' or '0'. A shape with increment '0' belongs to the same basis function of the closest previous shape with '1'.

`XL` to `ZU`: Lower and Upper x-, y-, and z-coordinates.

`Dir`: Surface normal direction.

`SDir`: Shape varying direction.

`Decay`: Positive sign for shape decaying toward the positive direction, and negative sign for decaying toward the negative direction.

`Ignore`: Reserved and currently not used in `caplet_solver`.

Examples can be found under folder `caplet_solver/example`.



####`.qui`

Files are originally used for FASTCAP. `caplet_geo` can also load and visualize them in 3D. The format is as the following:

    0
    shape conductor_label x1 y1 z1 x2 y2 z2 ...
    .
    .

The first line led by 0 is always ignored. From the second line to the last, each line defines a flat shape or a piecewise constant basis function. 

`shape`: `Q` for quadrilateral and `T` for triangular shapes.
`conductor_label`: arbitrary single word label. The common choice is integer indices.
`xn yn zn`: x-, y-, and z-coordinates for a single vertex. For `Q`, four such triplets follow `conductor_label`. For `T`, there are three triplets instead. The vertices should be recorded either clockwise or counter-clockwise.
In Manhattan geometries, we only use `Q` for piecewise constant basis functions.

More information can be found in [FASTCAP](http://www.rle.mit.edu/cpg/research_codes.htm) manual (`ug.tex` under `doc`).
