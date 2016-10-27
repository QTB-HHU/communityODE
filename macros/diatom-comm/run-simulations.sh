# Enter virtual environment
source ../../code/python/bin/activate
# Run community:
	# -j: with complete-media.json parameters
	# -f: in complete media 
	# -d: with diatom
	# -b: with pseudoalteromonas
	# -e: with flavobacteria
	# -a: with alteromonas
	# -c: with pseudomonas
	# -i: terminate simulation after 3000 integration steps
	# -l: terminate simulation at time 300
python communityODESys.py -j complete-media.json  -f -d -b -e -a -c -i 3000 -l 300
# Exit virtual environment
deactivate
