# PLONK in Plain English

Credits : [0xSachin](https://twitter.com/0xSachinK), [madab](https://twitter.com/iammadab), [mikefrancis](https://twitter.com/only1franchesco)

This is probably the simplest way to get to know how PLONK works. Just a few prereq :
- [ZKP Notes](https://github.com/0xSachinK/zkp-mooc-notes)
- [PLONK explainer](https://www.youtube.com/watch?v=A0oZVEXav24) by Dan Boneh 
- [Polynomial Commitments](https://www.youtube.com/watch?v=WyT5KkKBJUw) by Yupeng Zhang 

## Details to consider
- This implementation is not optimized for zero-knowledge. The parts of PlonK that are responsible for ensuring strong privacy are left out of this implementation.
- This implementation does not support Lookups & custom gates.


## Part 1 - The age of trusted setup

Ref : Plonk by Hand [Part 1](https://research.metastate.dev/plonk-by-hand-part-1/)

Let's start with our [setup.py](https://github.com/nullity00/plonkathon/blob/main/setup.py). 

### ``def from_file(cls, filename):``

- We start by reading a binary file which happens to be ``powersOfTau28_hez_final_11.ptau``. If you're curious to know how it looks, look at [ptau.txt](https://github.com/nullity00/plonkathon/blob/main/ptau.txt). If that didn't make sense to you, it's fine. We convert those binary vaues to integer values & then operate on them (look at [ptau.json](https://github.com/nullity00/plonkathon/blob/main/ptau.json)). 
- Byte 60 of the ptau file gives you the base-2 log of the value of ``powers``
-   ```
        contents[60] = log_2(d)
        2 ** contents[60] = d 
        d = powers ( in this case, d = 2^11 = 2048)
    ```
- According to the PLONK protocol paper, a circuit with ngates requires an SRS with at least n + 5 elements
- The elements in the ptau file have to be in the form of [G, xG, x2G, .... , H, xH, x2H, ...] where G & H are the generators of two subgroups G1 & G2 we intend to use in our pairing.
- We then extract G1 points which start at byte 80, by converting the bytes into integers. Make sure that none of the values are greater than the field modulus (prime p of a Field Fp) of the bn_128 curve using an assert ``assert max(values) < b.field_modulus``. 
- We know that ptau = [G, xG, x2G, .... , H, xH, x2H, ...] but looks like it is actually [yG, xyG, x2yG, .... , yH, xyH, x2yH, ...]. To extract the factor from the elements, we divide the first element by the generator ``yG/G`` (we know G, as the ptau is generated over bn128) & then remove the factor from all the elements by dividing it by y.
- We're done with extracting G1 points !
- Now, search for the starting point of G2 elements. We know that the first point in G2 elements is the generator H.
- ``ptau = [... (at 60)(log_2(d)), .... (at 80)(G), xG, .... (at 80 + 32 * powers * 2 - 1)(x^powers * G1) , ... , H, xH,.. G2 points ...]``
- The generator H (or G2) is 
```
G2 = (
    FQ2([
        10857046999023057135944570762232829481370756359578518086990519993285655852781,
        11559732032986387107991004021392285783925812861821192530917403151452391805634,
    ]),
    FQ2([
        8495653923123431417604973247489272438418190587263600148770280649306958101930,
        4082367875863433681332203403145435568316851327593401208105741076214120093531,
    ]),
)
```
- FQ2 is a class of field points in [py_ecc](https://github.com/ethereum/py_ecc/blob/master/py_ecc/fields/field_elements.py#L357)
- We then obtain the start of G2 points & return ``(powers_of_x, X2)`` where X2 is the starting point in G2.

### ``def commit(self, values: Polynomial) -> G1Point``

- This commit function requires the polynomials to be in Lagrange basis (This might not be the case always). 
- Run inverse FFT to convert values from Lagrange basis to monomial basis.
- Compute elliptic curve linear combination of setup with values ``[(G, a1), (xG, a2), (x2G, a3), ...]`` to return an elliptic curve point.

### ``def verification_key(self, pk: CommonPreprocessedInput) -> VerificationKey``

- To compute the verification key, we need the following :
```
class CommonPreprocessedInput:
    """Common preprocessed input"""

    group_order: int
    # q_M(X) multiplication selector polynomial
    QM: Polynomial
    # q_L(X) left selector polynomial
    QL: Polynomial
    # q_R(X) right selector polynomial
    QR: Polynomial
    # q_O(X) output selector polynomial
    QO: Polynomial
    # q_C(X) constants selector polynomial
    QC: Polynomial
    # S_σ1(X) first permutation polynomial S_σ1(X)
    S1: Polynomial
    # S_σ2(X) second permutation polynomial S_σ2(X)
    S2: Polynomial
    # S_σ3(X) third permutation polynomial S_σ3(X)
    S3: Polynomial
```
- We commit to all of the above polynomials & return the verification key.

## Part 2 - The never ending rounds

Ref : Plonk by Hand [Part 2](https://research.metastate.dev/plonk-by-hand-part-2-the-proof/)

Let's start with our [prover.py](https://github.com/nullity00/plonkathon/blob/main/prover.py). 

This is where we generate polynomials based on the inputs, outputs, selectors, constants. In theory, we are supposed to commit to these polynomials & the verifier queries these polynomials at a random point. In practice, to render the system non-interactive, we use Fiat-Shamir. We introduce random elements such as ``alpha``, ``beta``, ``gamma`` to evalute & then again commit to polynomials using random linear combination.

### Round 1

- Accumulate all the left, right, output witness values and compute `A`, `B`, `C` polynomials using lagrange interpolation.
- To generate a polynomial using values, we use the class `Polynomial` which takes in `n` witness values followed by `m - n` zeros (`m` is group order, `n` is total no. of wires) & the return basis of the polynomial
- Commit to these polynomials using the commit function in `setup.py`
- A full PLONK gate looks like ``Ql * A + Qr * B + Qo * C + Qm * A * C + Qc = PI`` where PI is the Witness Polynomial. We do a sanity check ``Ql * A + Qr * B + Qo * C + Qm * A * C + Qc - PI = Polynomial(0)`` to see if the polynomials were evaluated correctly.

### Round 2

- Round 2 is about committing to a single polynomial Z (permutation polynomial) which happens to encode all of the copy constraints. 
- We calculate `Z_values` by accumulating `N / D` where `N` is the product of RLCs of witness polynomials & the roots of unity & `D` is is the product of RLCs of witness polynomials & the selector polynomials
- The RLC function `self.rlc(val1, val2) = val_1 + self.beta * val_2 + gamma` where `beta` & `gamma` are random values we use to implement Fiat Shamir.
- Construct Z, Lagrange interpolation polynomial for Z_values & commit to Z

### Round 3

- Our goal is to compute the quotient polynomial `T(x)`, which will be of degree 3n + 5 for n gates. The polynomial `T(x)` encodes the majority of the information contained in our circuit and assignments all at once.
- 
```
    Open question:
    Why is the degree 3n + 5? Guess +5 is for ZK.
    If not zk, then 3 * (n + 1) - 1, for n gates
```
- We create cosets for polynomials using the `fft_expand` function. This function takes a random number, interpolates a polynomial at `[r, r * q, r * q**2 ... r * q**(4n-1)]` values where `r` is a random number . What is a coset btw ?
```
A coset, in the context of group theory, is a set that is formed by multiplying every element of a subgroup by a fixed element from the group.
```
- It is possible to find the quotient polynomial `T(x)` only if we satisfy these three conditions at the roots of unity.
    1. All gates are correct:
        `A * QL + B * QR + A * B * QM + C * QO + PI + QC = 0`
    2. The permutation accumulator is valid
        ```
        Z(wx) = Z(x) * N / D
        N = (rlc of A, X, 1) * (rlc of B, 2X, 1) * (rlc of C, 3X, 1) 
        D = (rlc of A, S1, 1) *(rlc of B, S2, 1) * (rlc of C, S3, 1)
        rlc = random linear combination: term_1 + beta * term2 + gamma * term3
        ```
    3. The permutation accumulator equals 1 at the start point
        ```
        (Z - 1) * L0 = 0
        L0 = Lagrange polynomial, equal at all roots of unity except 1
        Z = x^n - 1
        ```
- `T(x) = (Cosets of correct gates + cosets of permutation accumulator + Z_coset - 1 * L0 ) / Zh` 
- Since the degree of `T(x)` is too high, we split it up into `T1`, `T2`, `T3`
```
    T1 = Cosets of correct gates / Zh
    T2 = cosets of permutation accumulator / Zh
    T3 = Z_coset - 1 * L0 / Zh
```
- Do a sanity check to see if we've computed `T1`, `T2`, `T3` correctly using Barycentric evaluation. What is [Barycentric evaluation](https://hackmd.io/@vbuterin/barycentric_evaluation) ?
```
Given a polynomial expressed as a list of evaluations at roots of unity, evaluate it at x directly, without using an FFT to convert to coeffs first
```
- Commit to `T1`, `T2`, `T3`

### Round 4

- Compute A, B, C, S1, S2, Z polynomial at `zeta` (a random point we use to implement Fiat Shamir) using Barycentric evaluation.

### Round 5

- In Round 5, we create two large polynomials that combine all the polynomials we've been using so far and we output commitments to them.
- Compute the linearization polynomial, R & commit to it
- Compute the polynomials Wz & Wzw , commit to it

## PLONK in Rust

- [Plonk by Fingers](https://github.com/adria0/plonk-by-fingers)
- [Plonk by Hand](https://github.com/kcharbo3/Plonk-By-Hand)

## Now, try Goblin Plonk (This has nothing to do with PLONK)
![image](https://github.com/nullity00/plonkathon/assets/70228821/3ea7df5c-8d8b-44ef-9b63-75be31e04ee6)



