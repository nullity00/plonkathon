from utils import *
import py_ecc.bn128 as b
from curve import ec_lincomb, G1Point, G2Point
from compiler.program import CommonPreprocessedInput
from verifier import VerificationKey
from dataclasses import dataclass
from poly import Polynomial, Basis

# Recover the trusted setup from a file in the format used in
# https://github.com/iden3/snarkjs#7-prepare-phase-2
SETUP_FILE_G1_STARTPOS = 80
SETUP_FILE_POWERS_POS = 60


@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    # [x]₂ = xH, where H is a generator of G_2
    X2: G2Point

    @classmethod
    def from_file(cls, filename):
        contents = open(filename, "rb").read() 
        powers = 2 ** contents[SETUP_FILE_POWERS_POS] 

        # values = integer values of g1 points
        values = [
            # bytes -> integer in little endian
            int.from_bytes(contents[i : i + 32], "little") 

            # range(start, stop, step)
            # start = 80, stop = 80 + 32 * powers * 2, step = 32
            for i in range(
                SETUP_FILE_G1_STARTPOS, SETUP_FILE_G1_STARTPOS + 32 * powers * 2, 32
            )
        ]

        assert max(values) < b.field_modulus
        # The points are encoded in a weird encoding, where all x and y points
        # are multiplied by a factor (for montgomery optimization?). We can
        # extract the factor because we know the first point is the generator.
        factor = b.FQ(values[0]) / b.G1[0]
        values = [b.FQ(x) / factor for x in values]

        # range(powers) = [0, 1, 2, ..., powers - 1]
        # powers_of_x = [(values[0], values[1]), (values[2], values[3]), ...]
        # powers_of_x = [(G, x1G), (x2G, x3G), ...]
        powers_of_x = [(values[i * 2], values[i * 2 + 1]) for i in range(powers)]
        print("Extracted G1 side, X^1 point: {}".format(powers_of_x[1]))

        pos = SETUP_FILE_G1_STARTPOS + 32 * powers * 2

        # target = factor * 10857046999023057135944570762232829481370756359578518086990519993285655852781
        target = (factor * b.G2[0].coeffs[0]).n

        # checks where the target (G2 point) starts & assigns it to pos
        while pos < len(contents):
            v = int.from_bytes(contents[pos : pos + 32], "little")
            if v == target:
                break
            pos += 1
        print("Detected start of G2 side at byte {}".format(pos))
        # x2_encoding = contents[ pos + 128 : pos + 256 ]; len(x2_encoding) = 128
        X2_encoding = contents[pos + 32 * 4 : pos + 32 * 8]
        # X2_values = integer values of g2 points
        X2_values = [
            b.FQ(int.from_bytes(X2_encoding[i : i + 32], "little")) / factor
            for i in range(0, 128, 32)
        ]
        # x2 = first element of X2_values
        X2 = (b.FQ2(X2_values[:2]), b.FQ2(X2_values[2:]))
        assert b.is_on_curve(X2, b.b2)
        print("Extracted G2 side, X^1 point: {}".format(X2))

        # powers_of_x = [(G, x1G), (x2G, x3G), ...]
        # x2 = first element of X2_values (H)
        return cls(powers_of_x, X2)

    # Encodes the KZG commitment that evaluates to the given values in the group
    def commit(self, values: Polynomial) -> G1Point:
        assert values.basis == Basis.LAGRANGE

        # Run inverse FFT to convert values from Lagrange basis to monomial basis
        monomial_basis = values.ifft()
        # Optional: Check values size does not exceed maximum power setup can handle
        # monomial : 2x + 5x2
        assert len(monomial_basis.values) <= len(self.powers_of_x)
        # Compute linear combination of setup with values
        pairs = [] # pairs = [(G, a1), (xG, a2), (x2G, a3), ...]
        for i in range(len(monomial_basis.values)):
            pairs.append((self.powers_of_x[i], monomial_basis.values[i]))
        # computes G * a1 + xG * a2 + x2G * a3 + ...
        return ec_lincomb(pairs)

    # Generate the verification key for this program with the given setup
    def verification_key(self, pk: CommonPreprocessedInput) -> VerificationKey:
        # Create the appropriate VerificationKey object
        return VerificationKey(
            group_order=pk.group_order, 
            # commitment to multiplication selector polynomial
            Qm=self.commit(pk.QM), 
            # commitment to left selector polynomial
            Ql=self.commit(pk.QL), 
            # commitment to right selector polynomial
            Qr=self.commit(pk.QR), 
            # commitment to output selector polynomial
            Qo=self.commit(pk.QO), 
            # 	commitment to constants selector polynomial
            Qc=self.commit(pk.QC), 
            # commitment to the first permutation polynomial 
            S1=self.commit(pk.S1), 
            # commitment to the second permutation polynomial
            S2=self.commit(pk.S2), 
            # commitment to the third permutation polynomial
            S3=self.commit(pk.S3), 
            # X2 = xH, where H is a generator of G_2
            X_2=self.X2, 
            # nth root of unity, n - group order
            w=Scalar.root_of_unity(group_order=pk.group_order)
        )
