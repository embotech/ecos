#-------------------------------------------------------------------------------
# Makefile for all SuiteSparse packages
#-------------------------------------------------------------------------------

include SuiteSparse_config/SuiteSparse_config.mk

# Compile the default rules for each package
default:
	( cd SuiteSparse_config && $(MAKE) )
	( cd amd && $(MAKE) )
	( cd COLAMD && $(MAKE) library)
	( cd CHOLMOD && $(MAKE) library)

# Remove all files not in the original distribution
purge:
	- ( cd SuiteSparse_config && $(MAKE) purge )
	- ( cd amd && $(MAKE) purge )
	- ( cd COLAMD && $(MAKE) purge )
	- ( cd CHOLMOD && $(MAKE) purge )

# Remove all files not in the original distribution, but keep the libraries
clean:
	- ( cd SuiteSparse_config && $(MAKE) clean )
	- ( cd amd && $(MAKE) clean )
	- ( cd COLAMD && $(MAKE) clean )
	- ( cd CHOLMOD && $(MAKE) clean )

distclean: purge

