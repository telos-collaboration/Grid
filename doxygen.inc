# Copyright (C) 2004 Oren Ben-Kiki
# This file is distributed under the same terms as the Automake macro files.

# Generate automatic documentation using Doxygen. Goals and variables values
# are controlled by the various DX_COND_??? conditionals set by autoconf.
#
# The provided goals are:
# doxygen-doc: Generate all doxygen documentation.
# doxygen-run: Run doxygen, which will generate the documentation
#              (HTML, XML) but will not do the post
#              processing required for the rest of it (which is currently none).
#
# The following variable is intended for use in Makefile.am:
#
# DX_CLEANFILES = everything to clean.
#
# This is usually added to MOSTLYCLEANFILES.

## --------------------------------- ##
## Format-independent Doxygen rules. ##
## --------------------------------- ##

if DX_COND_doc

## ------------------------------- ##
## Rules specific for HTML output. ##
## ------------------------------- ##

if DX_COND_html

DX_CLEAN_HTML = @DX_DOCDIR@/html

endif DX_COND_html

## ------------------------------ ##
## Rules specific for XML output. ##
## ------------------------------ ##

if DX_COND_xml

DX_CLEAN_XML = @DX_DOCDIR@/xml

endif DX_COND_xml


.INTERMEDIATE: doxygen-run $(DX_PS_GOAL) $(DX_PDF_GOAL)

doxygen-run: @DX_DOCDIR@/@PACKAGE@.tag

doxygen-doc: doxygen-run

@DX_DOCDIR@/@PACKAGE@.tag: $(DX_CONFIG) $(pkginclude_HEADERS)
	rm -rf @DX_DOCDIR@
	$(DX_ENV) $(DX_DOXYGEN) $(srcdir)/$(DX_CONFIG)

DX_CLEANFILES = \
    @DX_DOCDIR@/@PACKAGE@.tag \
    -r \
    $(DX_CLEAN_HTML) \
    $(DX_CLEAN_XML)

endif DX_COND_doc
