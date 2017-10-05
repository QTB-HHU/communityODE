# Enter virtual environment
source ../../venv/bin/activate
# Run community:
	# -j: with complete or minimal media parameters
	# -f: in complete media 
	# -m: in minimal media 
	# -d: with diatom
	# -b: with pseudoalteromonas
	# -e: with flavobacteria
	# -a: with alteromonas
	# -c: with pseudomonas
	# -i: terminate simulation after 3000 integration steps
	# -l: terminate simulation at time 300
	# -P: plot paper Figure 5
python communityODESys.py -j complete-media-parset.json -f -d -b -e -a -c -i 3000 -l 300 -P
python communityODESys.py -j minimal-media-parset.json  -m -d -b -e -a -c -i 3000 -l 300 -P
# Exit virtual environment
deactivate
