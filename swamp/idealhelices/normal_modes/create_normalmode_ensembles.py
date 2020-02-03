import gemmi
import sys
import os


def clean_up(directory):
	for fname in os.listdir(directory):
		if "tmp_normalmode" in fname and os.path.isfile(os.path.join(directory, fname)):
			os.remove(os.path.join(directory, fname))

def run_nolb(pdbin):
	os.system('NOLB %s -o tmp_normalmode -a 0.2 -s 15' % pdbin)
	outfiles={}
	directory=os.getcwd()
	for fname in os.listdir(directory):
		if "tmp_normalmode" in fname and os.path.isfile(os.path.join(directory, fname)):
			mode = '%s_%s' % (fname.split("_")[2], fname.split("_")[3])
			outfiles[mode]=os.path.join(directory, fname)
	return outfiles


def get_models(pdbin, pdbout, models=(3,4,5,6,7,)):

	new_lines=[]
	count=0
	extract=False

	with open(pdbin, "r") as fhandle:
		for line in fhandle:
			if "MODEL" in line:
				if int(line.split()[1]) in models:
					extract=True
					count+=1
					new_lines.append("MODEL      %s\n" % count)
				else:
					extract=False
			elif extract:
				if "ENDMDL" in line:
					new_lines.append("ENDMDL\n")
				else:
					new_lines.append(line)
		new_lines.append("END\n")
	with open(pdbout, "w") as fhandle:
		for line in new_lines:
			fhandle.write(line)

	a = gemmi.read_structure(pdbout)
	for model in a:
		for chain in model:
			for residue in chain:
				residue.trim_to_alanine()
	a.write_minimal_pdb(pdbout)


os.chdir('/home/filo/PycharmProjects/CON-MOL/swamp/idealhelices/normal_modes')
outdir='/home/filo/PycharmProjects/CON-MOL/swamp/idealhelices/normal_modes'
polyalas={'polyala35': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala35.pdb', 'polyala10': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala10.pdb', 'polyala20': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala20.pdb', 'polyala15': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala15.pdb', 'polyala5': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala5.pdb', 'polyala40': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala40.pdb', 'polyala25': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala25.pdb', 'polyala30': '/home/filo/opt/CCP4/ccp4-7.0/share/ample/include/polyala30.pdb'}

for key in polyalas.keys():
	normal_modes = run_nolb(polyalas[key])
	for mode in normal_modes.keys():
		get_models(normal_modes[mode], os.path.join(outdir, '%s_%s' % (key, mode)))
	clean_up(os.getcwd())


