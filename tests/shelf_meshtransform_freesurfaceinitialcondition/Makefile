input: clean
	gmsh -format msh2 -2 src/square.geo -o square.msh
	../../tools/gmsh_mesh_transform.py 'shelflength = 500000, shelfslopeheight = 900000, minoceandepth = 100000, oceandepth = 1000000' 'x < shelflength' '(x, ((y + oceandepth)/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth) - oceandepth, z)' square.msh
	../../tools/gmsh_mesh_transform.py True '(x, 0.001 * y, z)' square.msh

clean:
	rm -f  *.ele *.edge *.node *.msh *.msh.bak *.vtu *.stat matrixdump* *.node.bak fluidity.*

run: input
	fluidity -v3 -l shelf.flml
