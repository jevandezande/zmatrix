'''
A module for converting between zmatrices and cartesian coordinates
'''

import math as m
import numpy as np

VERBOSE = True

class Converter:
	'''A coordinate converter class'''
	def __init__( self ):
		# Dictionary of the masses of elements indexed by element name; includes X for dummy atoms
		self.masses = { 'X': 0, 'Ac': 227.028, 'Al': 26.981539, 'Am': 243, 'Sb': 121.757, 'Ar': 39.948, 'As': 74.92159, 'At': 210, 'Ba': 137.327, 'Bk': 247, 'Be': 9.012182, 'Bi': 208.98037, 'Bh': 262, 'B': 10.811, 'Br': 79.904, 'Cd': 112.411, 'Ca': 40.078, 'Cf': 251, 'C': 12.011, 'Ce': 140.115, 'Cs': 132.90543, 'Cl': 35.4527, 'Cr': 51.9961, 'Co': 58.9332, 'Cu': 63.546, 'Cm': 247, 'Db': 262, 'Dy': 162.5, 'Es': 252, 'Er': 167.26, 'Eu': 151.965, 'Fm': 257, 'F': 18.9984032, 'Fr': 223, 'Gd': 157.25, 'Ga': 69.723, 'Ge': 72.61, 'Au': 196.96654, 'Hf': 178.49, 'Hs': 265, 'He': 4.002602, 'Ho': 164.93032, 'H': 1.00794, 'In': 114.82, 'I': 126.90447, 'Ir': 192.22, 'Fe': 55.847, 'Kr': 83.8, 'La': 138.9055, 'Lr': 262, 'Pb': 207.2, 'Li': 6.941, 'Lu': 174.967, 'Mg': 24.305, 'Mn': 54.93805, 'Mt': 266, 'Md': 258, 'Hg': 200.59, 'Mo': 95.94, 'Nd': 144.24, 'Ne': 20.1797, 'Np': 237.048, 'Ni': 58.6934, 'Nb': 92.90638, 'N': 14.00674, 'No': 259, 'Os': 190.2, 'O': 15.9994, 'Pd': 106.42, 'P': 30.973762, 'Pt': 195.08, 'Pu': 244, 'Po': 209, 'K': 39.0983, 'Pr': 140.90765, 'Pm': 145, 'Pa': 231.0359, 'Ra': 226.025, 'Rn': 222, 'Re': 186.207, 'Rh': 102.9055, 'Rb': 85.4678, 'Ru': 101.07, 'Rf': 261, 'Sm': 150.36, 'Sc': 44.95591, 'Sg': 263, 'Se': 78.96, 'Si': 28.0855, 'Ag': 107.8682, 'Na': 22.989768, 'Sr': 87.62, 'S': 32.066, 'Ta': 180.9479, 'Tc': 98, 'Te': 127.6, 'Tb': 158.92534, 'Tl': 204.3833, 'Th': 232.0381, 'Tm': 168.93421, 'Sn': 118.71, 'Ti': 47.88, 'W': 183.85, 'U': 238.0289, 'V': 50.9415, 'Xe': 131.29, 'Yb': 173.04, 'Y': 88.90585, 'Zn': 65.39, 'Zr': 91.224 }
		self.total_mass = 0
		self.cartesian = []
		self.zmatrix = []
		
	def read_zmatrix( self, input_file='zmatrix.dat' ):
		'''Read the input zmatrix file (assumes no errors and no variables)'''
		'''The zmatrix is a list with each element formatted as follows
		[ name, [[ atom1, distance ], [ atom2, angle ], [ atom3, dihedral ]], mass ]
		The first three atoms have blank lists for the undefined coordinates'''
		self.zmatrix = []
		with open( input_file, 'r' ) as f:
			f.readline()
			f.readline()
			name = f.readline().strip()
			self.zmatrix.append( [ name, [], self.masses[name] ] )
			name, atom1, distance = f.readline().split()[:3]
			self.zmatrix.append( [ name,
								   [ [int(atom1) - 1, float(distance)], [], [] ],
								   self.masses[name] ] )
			name, atom1, distance, atom2, angle = f.readline().split()[:5]
			self.zmatrix.append( [ name,
								   [[ int(atom1) - 1, float(distance) ],
								    [int(atom2) - 1, m.radians( float(angle) ) ], []],
								   self.masses[name] ] )
			for line in f.readlines():
				# Get the components of each line, dropping anything extra
				name, atom1, distance, atom2, angle, atom3, dihedral = line.split()[:7]
				# convert to a base 0 indexing system and use radians
				atom = [ name,
						[ [int(atom1) - 1, float(distance) ],
						  [int(atom2) - 1, m.radians( float(angle) ) ],
						  [int(atom3) - 1, m.radians( float(dihedral) ) ] ],
						self.masses[name] ]
				self.zmatrix.append( atom )
		
		return self.zmatrix
	
	def read_cartesian( self, input_file='cartesian.dat' ):
		'''Read the cartesian coordinates file (assumes no errors)'''
		'''The cartesian coordiantes consist of a list of atoms formatted as follows
		[ name, np.array( [ x, y, z ] ), mass ]
		'''
		self.cartesian = []
		with open( input_file, 'r' ) as f:
			#Throw away the first two lines
			f.readline()
			f.readline()
			for line in f.readlines():
				coords = line.split()
				name = coords[0]
				position = []
				for i in coords[1:4]:
					position.append( float(i) )
				self.cartesian.append( [name, np.array( position ), self.masses[name]] )

		return self.cartesian

	def rotation_matrix( self, axis, angle ):
		'''Euler-Rodrigues formula for rotation matrix'''
		# Normalize the axis
		axis = axis/np.sqrt(np.dot(axis,axis))
		a = np.cos( angle/2 )
		b,c,d = -axis*np.sin(angle/2)
		return np.array( [[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
						[2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
						[2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]] )

	def add_first_three_to_cartesian( self ):
		'''The first three atoms in the zmatrix need to be treated differently'''
		# First atom
		name, coords, mass = self.zmatrix[0]
		self.cartesian = [[ name, np.array( [0, 0, 0] ), mass ]]

		# Second atom
		name, coords, mass = self.zmatrix[1]
		distance = coords[0][1]
		self.cartesian.append( [ name, np.array( [distance, 0, 0] ), self.masses[name] ] )

		# Third atom
		name, coords, mass =  self.zmatrix[2]
		atom1, atom2 = coords[:2]
		atom1, distance = atom1
		atom2, angle = atom2
		q = np.array( self.cartesian[atom1][1] ) # position of atom 1
		r = np.array( self.cartesian[atom2][1] ) # position of atom 2

		# Vector pointing from q to r
		a = r - q

		# Vector of length distance pointing along the x-axis
		d = distance*a/np.sqrt(np.dot(a,a))

		# Rotate d by the angle around the z-axis
		d = np.dot( self.rotation_matrix( [0,0,1], angle ), d )

		# Add d to the position of q to get the new coordinates of the atom
		p = q + d
		atom = [ name, p, self.masses[name] ]
		self.cartesian.append( atom )

	def add_atom_to_cartesian( self, coords ):
		'''Find the cartesian coordinates of the atom'''
		name, coords, mass = coords
		atom1, distance = coords[0]
		atom2, angle = coords[1]
		atom3, dihedral = coords[2]
		
		q = self.cartesian[atom1][1] # atom 1
		r = self.cartesian[atom2][1] # atom 2
		s = self.cartesian[atom3][1] # atom 3
		
		# Vector pointing from q to r
		a = r - q
		# Vector pointing from s to r
		b = r - s

		# Vector of length distance pointing from q to r
		d = distance*a/np.sqrt(np.dot(a,a))

		# Vector normal to plane defined by q,r,s
		normal = np.cross( a, b )
		# Rotate d by the angle around the normal to the plane defined by q,r,s
		d = np.dot( self.rotation_matrix( normal, angle ), d )
		
		# Rotate d around a by the dihedral
		d = np.dot( self.rotation_matrix( a, dihedral ), d )

		# Add d to the position of q to get the new coordinates of the atom
		p = q + d
		atom = [ name, p, mass ]
		
		self.cartesian.append( atom )

	def zmatrix_to_cartesian( self ):
		'''Convert the zmartix to cartesian coordinates'''
		# Deal with first three line separately
		self.add_first_three_to_cartesian()

		for i in range( 3, len(self.zmatrix) ):
			self.add_atom_to_cartesian( self.zmatrix[i] )

		self.remove_dummy_atoms()

		self.center_cartesian()
		
		return self.cartesian
	
	def add_first_three_to_zmatrix( self ):
		'''The first three atoms need to be treated differently'''
		# First atom
		self.zmatrix = []
		name, position, mass = self.cartesian[0]
		self.zmatrix.append( [name, [[],[],[]], mass] )
		
		# Second atom
		name, position, mass = self.cartesian[1]
		atom1 = self.cartesian[0]
		pos1 = atom1[1]
		q = pos1 - position
		distance = m.sqrt( np.dot( q, q ) )
		self.zmatrix.append( [name, [[0,distance],[],[]], mass] )

		# Third atom
		name, position, mass = self.cartesian[2]
		atom1, atom2 = self.cartesian[:2]
		pos1, pos2 = atom1[1], atom2[1]
		q = pos1 - position
		r = pos2 - pos1
		q_u = q / np.sqrt( np.dot( q, q ) )
		r_u = r / np.sqrt( np.dot( r, r ) )
		distance = np.sqrt( np.dot( q, q ) )
		# Angle between a and b = acos( dot product of the unit vectors )
		angle = m.acos( np.dot( -q_u, r_u ) )
		self.zmatrix.append( [name, [ [ 0, distance ], [ 1, np.degrees( angle ) ], [] ], mass ] )
		
	def add_atom_to_zmatrix( self, i, line ):
		'''Generates an atom for the zmatrix
		(assumes that three previous atoms have been placed in the cartesian coordiantes)'''
		name, position, mass = line
		atom1, atom2, atom3 = self.cartesian[:3]
		pos1, pos2, pos3 = atom1[1], atom2[1], atom3[1]
		# Create vectors pointing from one atom to the next
		q = pos1 - position
		r = pos2 - pos1
		s = pos3 - pos2
		position_u = position / np.sqrt( np.dot( position, position ) )
		# Create unit vectors
		q_u = q / np.sqrt( np.dot( q, q ) )
		r_u = r / np.sqrt( np.dot( r, r ) )
		s_u = s / np.sqrt( np.dot( s, s ) )
		distance = np.sqrt( np.dot( q, q ) )
		# Angle between a and b = acos( dot( a, b ) / ( |a| |b| ) )
		angle = m.acos( np.dot( -q_u, r_u ) )
		angle_123 = m.acos( np.dot( -r_u, s_u ) )
		# Dihedral angle = acos( dot( normal_vec1, normal_vec2 ) / ( |normal_vec1| |normal_vec2| ) )
		plane1 = np.cross( q, r )
		plane2 = np.cross( r, s )
		dihedral = m.acos( np.dot( plane1, plane2 ) / ( np.sqrt( np.dot( plane1, plane1 ) ) * np.sqrt( np.dot( plane2, plane2 ) ) ) )
		# Convert to signed dihedral angle
		if np.dot( np.cross( plane1, plane2 ), r_u ) < 0:
			dihedral = -dihedral

		coords = [ [0, distance],[1, np.degrees( angle )],[2, np.degrees( dihedral )] ]
		atom = [ name, coords, mass ]
		self.zmatrix.append( atom )

	def cartesian_to_zmatrix( self ):
		'''Convert the cartesian coordinates to a zmatrix'''
		self.add_first_three_to_zmatrix()
		for i in range( 3, len(self.cartesian) ):
			line = self.cartesian[i]
			self.add_atom_to_zmatrix( i, line )

		return self.zmatrix

	def remove_dummy_atoms( self ):
		'''Delete any dummy atoms that may have been placed in the calculated cartesian coordinates'''
		new_cartesian = []
		for line in self.cartesian:
			if not line[0] == 'X':
				new_cartesian.append( line )
		self.cartesian = new_cartesian	

	def center_cartesian( self ):
		'''Find the center of mass and move it to the origin'''
		self.total_mass = 0.0
		center_of_mass = np.array( [ 0.0, 0.0, 0.0 ] )
		for atom in self.cartesian:
			mass = atom[2]
			self.total_mass += mass
			center_of_mass += atom[1]*mass
		center_of_mass = center_of_mass / float(self.total_mass)

		# Translate each atom by the center of mass
		for atom in self.cartesian:
			atom[1] = atom[1] - center_of_mass
	
	def cartesian_radians_to_degrees( self ):
		for atom in self.cartesian:
			atom[1][1][1] = np.degrees( atom[1][1][1] )  
			atom[1][2][1] = np.degrees( atom[1][2][1] )

	def output_cartesian( self, output_file='cartesian.dat' ):
		'''Output the cartesian coordinates of the file'''
		with open( output_file, 'w' ) as f:
			f.write( str( len(self.cartesian) ) )
			f.write( '\n\n' )
			for line in self.cartesian:
				name, position, mass = line 
				f.write( name + '\t' )
				f.write( '\t'.join( str(x) for x in position ) )
				f.write( '\n' )

	def print_cartesian( self ):
		'''Print the cartesian coordinates'''
		for line in self.cartesian:
			print line[0] + '\t' + '\t'.join( str(x) for x in line[1] )

	def output_zmatrix( self, output_file='zmatrix.dat' ):
		'''Output the zmatrix to the file'''
		with open( output_file, 'w' ) as f:
			f.write( '#ZMATRIX\n#\n' )
			for line in self.zmatrix:
				name, position, mass = line 
				f.write( name )
				for i in position:
					for j in range( 0, len(i), 2 ):
						f.write( '\t' + str(i[j]+1) + '\t' + str(i[j+1]) )
				f.write( '\n' )

	def print_zmatrix( self ):
		'''Print the zmatrix'''
		for line in self.zmatrix:
			print line[0] + '\t' + '\t'.join( str(x) for x in line[1] )

	def run_zmatrix( self, input_file='zmatrix.dat', output_file='cartesian.dat' ):
		'''Read in the zmatrix, converts it to cartesian, and outputs it to a file'''
		self.read_zmatrix( input_file )
		
		self.zmatrix_to_cartesian()

		self.output_cartesian( output_file )

		return 0
	
	def run_cartesian( self, input_file='cartesian.dat', output_file='zmatrix.dat' ):
		'''Read in the cartesian coordinates, convert to cartesian, and output the file'''
		self.read_cartesian( input_file )

		self.cartesian_to_zmatrix()

		self.output_zmatrix( output_file )

		return 0
