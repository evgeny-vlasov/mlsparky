# XSD validation
# 
# Author: Evgeny
###############################################################################



xsdLoc <- "http://nmrml.org/schema/v1.0.rc1/nmrML.xsd"
# validation
## internal R objects for validation
doc <- xmlInternalTreeParse(fileLoc)
tf <- xmlSchemaParse(xsdLoc)

## validation
valid <- xmlSchemaValidate(tf, doc)

valid 

