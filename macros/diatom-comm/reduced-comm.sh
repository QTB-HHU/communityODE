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
python communityODESys.py -j complete-media-parset.json -f -d -b -e -a -i 3000 -l 300 -P -t ' without P'
python communityODESys.py -j minimal-media-parset.json  -m -d -b -e -a -i 3000 -l 300 -P -t ' without P'
python communityODESys.py -j complete-media-parset.json -f -d -b -e -c -i 3000 -l 300 -P -t ' without A'
python communityODESys.py -j minimal-media-parset.json  -m -d -b -e -c -i 3000 -l 300 -P -t ' without A'
python communityODESys.py -j complete-media-parset.json -f -d -b -a -c -i 3000 -l 300 -P -t ' without F'
python communityODESys.py -j minimal-media-parset.json  -m -d -b -a -c -i 3000 -l 300 -P -t ' without F'
python communityODESys.py -j complete-media-parset.json -f -d -e -a -c -i 3000 -l 300 -P -t ' without PA'
python communityODESys.py -j minimal-media-parset.json  -m -d -e -a -c -i 3000 -l 300 -P -t ' without PA'
python communityODESys.py -j complete-media-parset.json -f -d -i 3000 -l 300 -P -t ', axenic'
python communityODESys.py -j minimal-media-parset.json  -m -d -i 3000 -l 300 -P -t ', axenic'
python communityODESys.py -j complete-media-parset.json -f -d -a -i 3000 -l 300 -P -t ', with A only'
python communityODESys.py -j minimal-media-parset.json  -m -d -a -i 3000 -l 300 -P -t ', with A only'
# Exit virtual environment
deactivate
