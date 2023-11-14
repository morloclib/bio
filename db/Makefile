all:
	morloc make main.loc

.PHONY: clean
clean:
	rm -rf nexus.py poo* __pycache__

.PHONY: test
test:
	./nexus.py searchEntrez '{"email":"zbwrnz@gmail.com","db":"nuccore","mindate":"2023/01/01","maxdate":"2023/01/31","retmax":100}' '"Influenza+A+Virus[Organism]+H3N2"'

.PHONY: testrec
testrec:
	./nexus.py nucleotideAccessionToJson '{"email":"zbwrnz@gmail.com"}' '["OQ335855.1","ON965379.1","ON965378.1"]'
