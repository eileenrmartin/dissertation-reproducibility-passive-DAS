include Makefile-ER

ER:
	make all-ER-figs

# can be run on ***which sep machine*****
CR-SEP:
	make all-CR-SEP-figs

# can be run on cees-mazama, cees-tool-7, cees-tool-8.stanford.edu
CR-CEES:
	make all-CR-CEES-figs

# can be run on snowbear.lbl.gov (don't forget regular password followed by Google authenticator token)
CR-Bears:
	make all-CR-Bears-figs