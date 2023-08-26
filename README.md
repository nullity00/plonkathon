# PlonKathon - PLONK in Plain English

Credits : [0xSachin](https://twitter.com/0xSachinK), [madab](https://twitter.com/iammadab), [mikefrancis](https://twitter.com/only1franchesco)

This is probably the simplest way to get to know how PLONK works. Just a few prereq : 
- [PLONK explainer](https://www.youtube.com/watch?v=A0oZVEXav24) by Dan Boneh | [Notes](https://hackmd.io/mNr5tcxtRTiStp9DC5oSXg)
- [Polynomial Commitments](https://www.youtube.com/watch?v=WyT5KkKBJUw) | [Notes](https://hackmd.io/S9L9JGWUQ2W-2-NA24-5KQ)


## PLONK - The age of trusted setup

Ref : Plonk by Hand [Part 1](https://research.metastate.dev/plonk-by-hand-part-1/)

Let's start with our [setup.py](https://github.com/nullity00/plonkathon/blob/main/setup.py). 

### def from_file(cls, filename):

- We start by reading a binary file which happens to be ``powersOfTau28_hez_final_11.ptau``. If you're curious to know how it looks, look at [ptau.txt](). If that didn't make sense to you, it's fine. We convert those binary vaues to integer values & then operate on them (look at [ptau.json]()). 
- Byte 60 of the ptau file gives you the base-2 log of the value of ``powers``
-   ```
        contents[60] = log_2(d)
        2 ** contents[60] = d 
        d = powers ( for this case, d = 2^11 = 2048)
    ```
- According to the PLONK protocol paper, a circuit with ngates requires an SRS with at least n + 5 elements
- The elements in the ptau file have to be in the form of [G, xG, x2G, .... , H, xH, x2H, ...] where G & H are the generators of two subgroups G1 & G2 we intend to use in our pairing.
- We then extract G1 points which start at byte 80, by converting the bytes into integers. Make sure that none of the values are greater than the field modulus (prime p of a Field Fp) of the bn_128 curve using an assert ``assert max(values) < b.field_modulus``. 
- We know that ptau = [G, xG, x2G, .... , H, xH, x2H, ...] but looks like it is actually [yG, xyG, x2yG, .... , yH, xyH, x2yH, ...]. To extract the factor from the elements, we divide the first element by the generator ``yG/G`` (we know G, as the ptau is generated over bn128) & then remove the factor from all the elements by dividing it by y.
- We're done with extracting G1 points !
- Now, search for the starting point of G2 elements. We know that the first point in G2 elements is the generator H.
- ``ptau = [... (at 60)(log_2(d)), .... (at 80)(G), xG, .... (at 80 + 32 * powers * 2 - 1)(x^powers * G1) , ... , H, xH,.. G2 points ...]``
- The generator H (or G2) is 
```G2 = (
    FQ2([
        10857046999023057135944570762232829481370756359578518086990519993285655852781,
        11559732032986387107991004021392285783925812861821192530917403151452391805634,
    ]),
    FQ2([
        8495653923123431417604973247489272438418190587263600148770280649306958101930,
        4082367875863433681332203403145435568316851327593401208105741076214120093531,
    ]),
)```
- FQ2 is a class of field points in ``[py_ecc](https://github.com/ethereum/py_ecc/blob/master/py_ecc/fields/field_elements.py#L357)``
- We then obtain the start of G2 points & return ``(powers_of_x, X2)`` where X2 is the starting point in G2.

### def commit(self, values: Polynomial) -> G1Point:

- 

