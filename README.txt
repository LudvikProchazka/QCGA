  Please note that this is a school project created primarily to learn C++ and to perform calculations in geometric algebra. The goal was to make the calculations feel more human-readable
and as easy to debug as possible, even at the expense of speed and performance. As a result, the code architecture is poor — but it works. :)
With the knowledge I have gained over the years, I would now choose a different design. However, since implementing high-dimensional geometric algebras is a challenge in itself, I have decided to leave it as it is for now.
  Calculations are done by directly writing into the main function. There are few examples prepared, how to use it. In QCGA.h, there are defines for 1, e1, e2, ... , e15, as well as for the null basis ei1,eo1,...ei6,eo6.
Those defines are the building blocks for all multivectors to be used in the main function. Please use some functions carefully, as they may cause exceptions if they are used with not suitible multivector (because of bad design).
Thus, a user has to know exactly what he is doing.

What to be aware of:
-QCGA object having scalar part has to be of the form "scalar" * one + (rest of the multivector), since "one" is an object representing 1 as an algebra's identity.
    Eg: QCGA vec = 5 * one + e1 * e2;
-The operator precedence is sometimes wild so when something doesn't work you may have not enough parentheses
-Examples in Example.h serve to some extend as unit tests.
-Both the design and implementation may be improved, because at the time the code was written, I was very new to C++ and quite new to programming in general.

Have fun and good luck exploring it :) If you have any idea or question, feel free to reach me: 209458@vutbr.cz
