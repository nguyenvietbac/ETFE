import time
import math
import random as rd
import numpy as np
from numpy import linalg as LA


start_time = time.time()
# initial data
file_out_name = 'semi.dat'
density = 0.05
c_length =  1000
n_chains = 10   
b_length = 1.53
seed = 9062
mass = 14.02

sqrt_2_2 = np.sqrt(2)/2
n_oc_sites = c_length*n_chains
t_sites = n_oc_sites/density
n_incs = int(np.cbrt(t_sites/4)) + 1



def find_empty_sites(curr_site, sites):
	# d = np.zeros((3, 3, 3, 4), dtype=np.int)
	dx = np.zeros(3, dtype=np.int)
	dy = np.zeros(3, dtype=np.int)
	dz = np.zeros(3, dtype=np.int)
	# m = np.zeros(4, dtype=np.int)
	a = np.zeros(3)
	aa = np.zeros(3)
	va = np.zeros(3)
	empty_sites_around = []

	m = curr_site['x3']
	dx[1] = curr_site['x0']
	dy[1] = curr_site['x1']
	dz[1] = curr_site['x2']
	dx[2] = curr_site['x0']+1
	dy[2] = curr_site['x1']+1
	dz[2] = curr_site['x2']+1
	dx[0] = curr_site['x0']-1
	dy[0] = curr_site['x1']-1
	dz[0] = curr_site['x2']-1

	if dx[2]==n_incs:
		dx[2] = 0
	if dy[2]==n_incs:
		dy[2] = 0
	if dz[2]==n_incs:
		dz[2] = 0
	if dx[0] < 0:
		dx[0] = n_incs - 1 
	if dy[0] < 0:
		dy[0] = n_incs - 1
	if dz[0] < 0:
		dz[0] = n_incs - 1

	if m==0:
		a[0] = 0
		a[1] = 0
		a[2] = 0
	elif m==1:
		a[0] = 0.5
		a[1] = 0.5
		a[2] = 0
	elif m==2:
		a[0] = 0.5
		a[1] = 0
		a[2] = 0.5
	elif m==3:
		a[0] = 0
		a[1] = 0.5
		a[2] = 0.5

	for mm in range(4):
		for x in range(3):
			for y in range(3):
				for z in range(3):
					if mm==0:
						aa[0] = 0
						aa[1] = 0
						aa[2] = 0
					elif mm==1:
						aa[0] = 0.5
						aa[1] = 0.5
						aa[2] = 0
					elif mm==2:
						aa[0] = 0.5
						aa[1] = 0
						aa[2] = 0.5
					elif mm==3:
						aa[0] = 0
						aa[1] = 0.5
						aa[2] = 0.5
					aa[0] += x - 1	
					aa[1] += y - 1
					aa[2] += z - 1
					va[0] = aa[0] - a[0]
					va[1] = aa[1] - a[1]
					va[2] = aa[2] - a[2]
					# print(aa)
					va_mag = LA.norm(va)
					if 0.5<va_mag<=sqrt_2_2+0.01:
						if sites[dx[x]][dy[y]][dz[z]][mm] == 0:
							temp = {}
							temp['x0'] = dx[x]
							temp['x1'] = dy[y]
							temp['x2'] = dz[z]
							temp['x3'] = mm
							temp['chain'] = curr_site['chain']
							temp['type'] = curr_site['type']
							temp['vec'] = np.copy(va)
							empty_sites_around.append(temp)

	return empty_sites_around

def main():
	# initial lacttic
	sites = np.zeros((n_incs,n_incs,n_incs,4), dtype=np.int)
	past_sites = np.copy(sites)
	future_sites = np.zeros((n_incs,n_incs,n_incs,4), dtype=np.int)
	# initial atoms list in lacttic
	atoms_list = []
	x = 0
	while x < n_oc_sites:
		temp = {}
		temp['x0'] = rd.randint(0,n_incs - 1)
		temp['x1'] = rd.randint(0,n_incs - 1)
		temp['x2'] = rd.randint(0,n_incs - 1)
		temp['x3'] = rd.randint(0,3)
		temp['vec'] = []
		for y in atoms_list:
			if temp == y:
				temp = 1
				break 
		if temp == 1:
			continue

		sites[temp['x0']][temp['x1']][temp['x2']][temp['x3']] = 1
		temp['type'] = 1
		atoms_list.append(temp)
		x = x + 1

		# print(type(temp))
	# print(sites)
	# create the chains
	x = 0
	ini_chain = 1
	chains_num = 1
	length_chain = 0
	list_atoms_in_chains = []
	while x < n_oc_sites:
		rd_temp = rd.randint(0,len(atoms_list) - 1)
		found = 0
		if ini_chain == 1:
			# at_temp = {}
			atoms_list[rd_temp]['chain'] = chains_num
			atoms_list[rd_temp]['vec'] = [0, 0, 0]
			at_temp = dict(atoms_list[rd_temp])
			first_site = atoms_list[rd_temp]
			del atoms_list[rd_temp]
			rd_temp = rd.randint(0,len(atoms_list) - 1)
			empty_sites_around = find_empty_sites(first_site, sites)
			while found != 1:
				if empty_sites_around == []:
					atoms_list.append(at_temp)
					ini_chain = 1
					break

				temp1 = rd.randint(0, len(empty_sites_around) - 1)
				second_site	= empty_sites_around[temp1]
				del empty_sites_around[temp1]
				va1 = second_site['vec']

				empty_sites_around2 = find_empty_sites(second_site, sites)
				for check_site in empty_sites_around2:
					vb1 = check_site['vec']
					angle = math.acos((np.dot(va1, vb1))/(LA.norm(va1)*LA.norm(vb1)))
					angle = round(angle*180/math.pi)
					if 40<angle<70:
						found = 1
						ini_chain = 0
						list_atoms_in_chains.append(at_temp)
						x += 1
						list_atoms_in_chains.append(second_site)
						x += 1
						sites[atoms_list[rd_temp]['x0']][atoms_list[rd_temp]['x1']][atoms_list[rd_temp]['x2']][atoms_list[rd_temp]['x3']] = 0
						sites[second_site['x0']][second_site['x1']][second_site['x2']][second_site['x3']] = 1
						length_chain = length_chain + 2
						break

		elif ini_chain == 0:
			empty_sites_around = find_empty_sites(list_atoms_in_chains[-1], sites)
			# print(empty_sites_around)
			while found != 1:
				if empty_sites_around == []:
					length_chain = 0
					chains_num += 1
					ini_chain = 1
					break
				temp1 = rd.randint(0, len(empty_sites_around) - 1)
				next_site = empty_sites_around[temp1]
				del empty_sites_around[temp1]
				va1 = next_site['vec']
				empty_sites_around2 = find_empty_sites(next_site, sites)
				for check_site in empty_sites_around2:
					vb1 = check_site['vec']
					angle = math.acos((np.dot(va1, vb1))/(LA.norm(va1)*LA.norm(vb1)))
					angle = round(angle*180/math.pi)
					# print(angle)
					if 40<angle<70:
						found = 1
						list_atoms_in_chains.append(next_site)
						x += 1
						length_chain += 1
						sites[atoms_list[rd_temp]['x0']][atoms_list[rd_temp]['x1']][atoms_list[rd_temp]['x2']][atoms_list[rd_temp]['x3']] = 0
						sites[next_site['x0']][next_site['x1']][next_site['x2']][next_site['x3']] = 1						
						break

			# break

						
		if length_chain == c_length:
			length_chain = 0
			chains_num += 1
			ini_chain = 1
	count = 0 
	for x in range(len(list_atoms_in_chains)):
		for y in range(x+1,len(list_atoms_in_chains)):
			if list_atoms_in_chains[x]['x0'] == list_atoms_in_chains[y]['x0']:
				if list_atoms_in_chains[x]['x1'] == list_atoms_in_chains[y]['x1']:
					if list_atoms_in_chains[x]['x2'] == list_atoms_in_chains[y]['x2']:
						if list_atoms_in_chains[x]['x3'] == list_atoms_in_chains[y]['x3']:
							count = count + 1
							print(list_atoms_in_chains[x])
	# for y in list_atoms_in_chains:
	# 	# print(y)
	# 	print(type(y))
	# # # 	# break

	chains_num = chains_num - 1
	# write in file

	file = open("PE_python.dat", "w" )
	file.write("# Model for PE by python copyright NguyenVietBac\n")
	file.write("\n")
	file.write("%10i     atoms\n" % (n_oc_sites))
	file.write("%10i     bonds\n" % (n_oc_sites - chains_num))
	file.write("%10i     angles\n" % (n_oc_sites - 2*chains_num))
	file.write("%10i     dihedrals\n" % (n_oc_sites - 3*chains_num))
	file.write("\n")
	file.write("%10i     atom types\n" % 1)
	file.write("%10i     bond types\n" % 1)
	file.write("%10i     angle types\n" % 1)
	file.write("%10i     dihedral types\n" % 1)
	file.write("\n")
	file.write("%10.4f%10.4f xlo xhi\n" % (0.0, n_incs*b_length*np.sqrt(2)))
	file.write("%10.4f%10.4f ylo yhi\n" % (0.0, n_incs*b_length*np.sqrt(2)))
	file.write("%10.4f%10.4f zlo zhi\n" % (0.0, n_incs*b_length*np.sqrt(2)))
	file.write("\n")
	file.write("Masses\n")
	file.write("\n")
	file.write("%10i %14.2f\n" % (1, mass))
	file.write("\n")
	file.write("Atoms\n")
	file.write("\n")
	a = np.zeros(3)
	count = 0
	for x in list_atoms_in_chains:
		xx = dict(x)
		# print(type(x))
		if type(x) == 'numpy.ndarray':
			print(numpy.ndarray)
		if xx['x3'] == 0:
			a[0] = 0
			a[1] = 0
			a[2] = 0
		elif xx['x3'] == 1:
			a[0] = 0.5
			a[1] = 0.5
			a[2] = 0
		elif xx['x3'] == 2:
			a[0] = 0.5
			a[1] = 0
			a[2] = 0.5
		elif xx['x3'] == 3:
			a[0] = 0
			a[1] = 0.5
			a[2] = 0.5

		xx['x0'] = b_length*np.sqrt(2)*(xx['x0'] + a[0])
		xx['x1'] = b_length*np.sqrt(2)*(xx['x1'] + a[1])
		xx['x2'] = b_length*np.sqrt(2)*(xx['x2'] + a[2])
		count += 1 
		file.write("%10i%10i%10i%10.4f%10.4f%10.4f\n" % (count, xx['chain'], xx['type'], xx['x0'], xx['x1'], xx['x2']))

	file.write("\n")
	file.write("Bonds \n")
	file.write("\n")
	count = 1
	atom_1 = 1
	atom_2 = 2
	for x in range(n_chains):
		for y in range(c_length - 1):
			file.write("%10i%10i%10i%10i\n" % (count, 1, atom_1, atom_2))
			count += 1
			atom_1 += 1
			atom_2 += 1
		atom_1 += 1
		atom_2 += 1
	file.write("\n")
	file.write("Angles \n")
	file.write("\n")
	count = 1
	atom_1 = 1
	atom_2 = 2
	atom_3 = 3
	for x in range(n_chains):
		for y in range(c_length - 2):
			file.write("%10i%10i%10i%10i%10i\n" % (count, 1, atom_1, atom_2, atom_3))
			count += 1
			atom_1 += 1
			atom_2 += 1 
			atom_3 += 1
		atom_1 += 2
		atom_2 += 2
		atom_3 += 2

	file.write("\n")
	file.write("Dihedrals \n")
	file.write("\n")
	count = 1
	atom_1 = 1
	atom_2 = 2
	atom_3 = 3
	atom_4 = 4

	for x in range(n_chains):
		for y in range(c_length - 3):
			file.write("%10i%10i%10i%10i%10i%10i\n" % (count, 1, atom_1, atom_2, atom_3, atom_4))
			count += 1
			atom_1 += 1
			atom_2 += 1
			atom_3 += 1
			atom_4 += 1
		atom_1 += 3
		atom_2 += 3
		atom_3 += 3
		atom_4 += 3


if __name__== "__main__":
	main()

end_time = time.time()

print(end_time - start_time)