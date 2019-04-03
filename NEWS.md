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
