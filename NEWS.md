<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.9 -- Release notes

### Implementation changes

	- Adapt collections to new TKCat

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.8 -- Release notes

### DESCRIPTION changes

   - Author information

### Implementation changes
   
   - Remove trailing spaces.
   
<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.7 -- Release notes

### Implementation changes
   
   - Adding "Likely pathogenic, Affects" to the *ClinVar_clinSigOrder*
   custom data file.

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.6 -- Release notes

### Data model changes
   
   - Fixed foreign key between varEntrez and entrezNames
   
### Implementation changes
   
   - Use of json data model

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.5 -- Release notes

### Data model changes
   
   - display_stop field of the ClinVar_varSeqLoc table is now **nullable**

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.4 -- Release notes

### Data model changes
   
   - Add collections

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.3 -- Release notes

### Data model changes

   - The sampleDescription field has been removed from the ClinVar_rcvaObservedIn
   table because it is not found in the xml source file anymore.
   - The assertion field of the ClinVar_ReferenceClinVarAssertion table
   has been set to NOT NULL
   - The type of the current field of the ClinVar_sourceFiles table has been
   set to Date.

### Implementation changes
   
   - Use of the "here" package
   - Write last update

<!----------------------------------------------------------------------------->
<!----------------------------------------------------------------------------->
## Version 0.2 -- Release notes

### Implementation changes

	- Tables are exported with " quotes for delimiting fields and escaped when needed by doubling them
