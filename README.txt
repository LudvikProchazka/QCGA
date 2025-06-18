This school project aims to help me with calculations in geometric algebra for quadratic surfaces (personal use). It is an extension of the conformal geometric algebra (CGA) and geometric algebra for conics (GAC) projects.
The aim is to manipulate multivectors in possibly intuitive way, not to be quick since with the dimension 2^15 it gets messy. 
What to be aware of:
-QCGA object having scalar part has to be of the form \scalar*one + (rest of the multivector)\ since "one" is an object representing 1 as an algebra's identity.
-The operator precedence is sometimes wild so when something doesn't work you may have not enough parentheses
-Examples in Example.h serve to some extend as unit tests.
-I had no idea someone will be interested in this project when I started :)
-Many improvements may be done - for example using std::map<int, double> where int will be 1-15 generating basis vectors + 0.
