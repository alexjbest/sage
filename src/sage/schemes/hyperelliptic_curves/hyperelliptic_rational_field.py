"""
Hyperelliptic curves over the rationals
"""
#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.abc

from sage.rings.padics.all import pAdicField

from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field
from sage.libs.pari.all import pari
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from . import hyperelliptic_generic
from .hyperelliptic_g2 import HyperellipticCurve_g2
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])


class HyperellipticCurve_rational_field(hyperelliptic_generic.HyperellipticCurve_generic,
                                        ProjectivePlaneCurve_field):

    def matrix_of_frobenius(self, p, prec=20):

        # BUG: should get this method from HyperellipticCurve_generic
        def my_chage_ring(self, R):
            from .constructor import HyperellipticCurve
            f, h = self._hyperelliptic_polynomials
            y = self._printing_ring.gen()
            x = self._printing_ring.base_ring().gen()
            return HyperellipticCurve(f.change_ring(R), h, "%s,%s"%(x,y))

        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        if isinstance(p, (sage.rings.abc.pAdicField, sage.rings.abc.pAdicRing)):
            K = p
        else:
            K = pAdicField(p, prec)
        frob_p, forms = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(my_chage_ring(self, K))
        return frob_p

    @cached_method
    def minimal_discriminant(self):
        r"""
        Return the minimal discriminant of this hyperelliptic curve as an integer.

        This method is cached.

        EXAMPLES::

            sage: F,G = EllipticCurve([0,0,1,-1,0]).hyperelliptic_polynomials()
            sage: E = HyperellipticCurve(F,G)
            sage: E.minimal_discriminant()
            37

            sage: F,G = EllipticCurve([0, -1, 1, -10, -20]).hyperelliptic_polynomials()
            sage: E = HyperellipticCurve(F,G)
            sage: E.minimal_discriminant()
            -161051

        LMFDB curve 113236.a.452944.1::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([0, 2, 12, 24, 16, 1]), R([0, 1, 1]))
            sage: C.minimal_discriminant()
            -452944
            sage: X = HyperellipticCurve(R([0, 8, 49, 98, 65, 4]))
            sage: X.minimal_discriminant()
            -452944

        LMFDB curve 3125.a.3125.1::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([256, 0, 0, 0, 0, 1]), R([]))
            sage: C.minimal_discriminant()
            3125
            sage: X = HyperellipticCurve(R([0,0,0,0,0,1]), R(1))
            sage: X.minimal_discriminant()
            3125

        LMFDB curve 50000.a.200000.1::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([1/256, 0, 0, 0, 0, 1]), R([]))
            sage: C.minimal_discriminant()
            200000
            sage: X = HyperellipticCurve(R([0,0,0,0,0,2]), R(1))
            sage: X.minimal_discriminant()
            200000

        LMFDB curve 400.a.409600.1::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([1, 0, 4, 0, 4, 0, 1]), R([]))
            sage: C.minimal_discriminant()
            -409600
            sage: X = HyperellipticCurve(R([1, 0, 4, 0, 4, 0, 1]))
            sage: X.minimal_discriminant()
            -409600

        Qing Lui's example from https://www.math.u-bordeaux.fr/~qliu/articles/algo_Weierstrass_mini-arxiv.pdf::

            sage: R.<x> = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(5^6*17^3*x^5, 2^4*11*13)
            sage: C.minimal_discriminant().factor()
            2^12 * 5^11 * 11^8 * 13^8 * 17^8

        A genus 3 example from https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.12604
        cross-referencing with the optional package
        https://alexjbest.github.io/cluster-pictures/cluster_pictures.html#sage_cluster_pictures.cluster_pictures.BYTree.minimal_discriminant::

            sage: C = HyperellipticCurve((x^3 - 7^15)*(x^2-7^6)*(x^3-7^3))
            sage: C.minimal_discriminant()
            590561633821394821581347105694543392818159361592405117014533630178101753364407699502733716603631767090783846400000000000000000000
            sage: C.minimal_discriminant().factor()
            2^64 * 3^26 * 5^20 * 7^24 * 13^10 * 19^10 * 43^10 * 181^10
            sage: (7,24) in C.minimal_discriminant().factor()
            True
            sage: (7,22) in HyperellipticCurve(7*(x^2+1)*(x^2+36)*(x^2+64)).minimal_discriminant().factor()

        An infinite family of examples::

            sage: p = 11
            sage: (p, 27) in HyperellipticCurve(p*(x^2-p^5)*(x^3-p^3)*((x-1)^3-p^9)).minimal_discriminant().factor()
            True
            sage: p = 13
            sage: (p, 27) in HyperellipticCurve(p*(x^2-p^5)*(x^3-p^3)*((x-1)^3-p^9)).minimal_discriminant().factor()
            True

        Curiously this holds even for even a couple of primes smaller than what is proved in the User's guide::

            sage: p = 5
            sage: (p, 27) in HyperellipticCurve(p*(x^2-p^5)*(x^3-p^3)*((x-1)^3-p^9)).minimal_discriminant().factor()
            True
            sage: p = 7
            sage: (p, 27) in HyperellipticCurve(p*(x^2-p^5)*(x^3-p^3)*((x-1)^3-p^9)).minimal_discriminant().factor()
            True


        """
        return pari(self.normalize_defining_polynomials().hyperelliptic_polynomials()
                   ).hyperellminimaldisc().sage()

    def rational_points(self, bound=None):
        r"""
        Find rational points on the hyperelliptic curve, all arguments are passed
        on to :meth:`sage.schemes.generic.algebraic_scheme.rational_points`.

        EXAMPLES:

        For the LMFDB genus 2 curve `932.a.3728.1 <https://www.lmfdb.org/Genus2Curve/Q/932/a/3728/1>`::

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, 1, 0, 1, -2, 1]), R([1]));
            sage: C.rational_points(bound=8)
            [(-1 : 2 : 1),
             (-1 : -3 : 1),
             (0 : 0 : 1),
             (0 : -1 : 1),
             (1 : 0 : 1),
             (1 : -1 : 1),
             (1/2 : -3/8 : 1),
             (1/2 : -5/8 : 1),
             (0 : 1 : 0)]

        Check that :trac:`29509` is fixed for the LMFDB genus 2 curve
        `169.a.169.1 <https://www.lmfdb.org/Genus2Curve/Q/169/a/169/1>`::

            sage: C = HyperellipticCurve(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]));
            sage: len(C.rational_points(bound=10))
            5
        
        ::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([1/256, 0, 0, 0, 0, 1]), R([]))
            sage: len(C.rational_points(bound=10))
            5

        """
        if not bound:
            raise TypeError("Unable to enumerate points over Q, please specify a height bound.")
        pariout = pari(self.hyperelliptic_polynomials()).hyperellratpoints(pari(bound)).sage()
        return [self(P) for P in pariout] + [self(0,1,0)]

    def is_locally_solvable(self, p) -> bool:
        r"""
        Whether the projective curve given by the affine model is solvable at prime `p`

        EXAMPLES::

            sage: R = PolynomialRing(QQ, "x")
            sage: C = HyperellipticCurve(R([1/256, 0, 0, 0, 0, 1]), R([]))
            sage: C.is_locally_solvable(5)
            True
            sage: C.is_locally_solvable(-1)
            Traceback (most recent call last):
            ...
            ValueError: curve must be of genus 2
        
        Lind and Reichardt's example of a curve everywhere locally solvable but with no rational points::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(2*x^4 - 34)
            sage: C.is_locally_solvable(17)
            True
            sage: C.is_locally_solvable(2)
            True

        ::

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([-56, 0, -75, 0, 15, 0, -1]), R([0, 1, 0, 1]));
            sage: C.is_locally_solvable(2)
            False
            sage: C.is_locally_solvable(3)
            True
            sage: C.is_locally_solvable(7)
            True

        """
        return pari(self.simplified_model().normalize_defining_polynomials().hyperelliptic_polynomials()[0]
                      ).hyperell_locally_solvable(ZZ(p).pari()).sage()

    def bad_primes(self):
        return

    def is_everywhere_locally_solvable(self, return_prime):
        r"""
        Whether the projective curve given by the affine model is solvable at prime `p`

        Warning: the Sage convention is that `rational_points` of a hyperelliptic curve
        are those of the singular model, thus it is possible that a hyperelliptic 

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(R([-56, 0, -75, 0, 15, 0, -1]), R([0, 1, 0, 1]))
            sage: C.is_everywhere_locally_solvable()
            False

        Lind and Reichardt's example of a curve everywhere locally solvable but with no rational points::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(2*x^4 - 34)
            sage: C.is_everywhere_locally_solvable()
            True

        """


    def reduction_type(self):
        r"""
        Compute the reduction type of a genus 2 hyperelliptic curve over `\QQ` using the
        ``genus2reduction`` program.

        EXAMPLES:

        X_1(13)::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]));
            sage: C.reduction_type()
            Reduction data about this proper smooth genus 2 curve:
                y^2 + (x^3 + x + 1)*y = x^5 + x^4
            A Minimal Equation:
                y^2 + (x^3 + x + 1)y = x^5 + x^4
            Minimal Discriminant: -169
            Conductor: 169
            Local Data:
                    p=13
                    (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                    reduction at p: [I{0}-II-0] page 159, (1), f=2

        Can't do higher genus yet::

            sage: x = polygen(QQ, 'x')
            sage: C = HyperellipticCurve(x^2+x, x^7+x^2+1)
            sage: C.reduction_type()
            Traceback (most recent call last):
            ...
            ValueError: curve must be of genus 2
        
        """
        if not isinstance(self, HyperellipticCurve_g2):
            raise ValueError('curve must be of genus 2')
        P,Q = self.hyperelliptic_polynomials()
        return genus2reduction(Q, P)

    def conductor(self, odd_part_only=False):
        r"""
        Compute the conductor of a genus 2 hyperelliptic curve over `\QQ`.

        EXAMPLES:

        X_1(13)::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]));
            sage: C.conductor()
            169
            sage: X = HyperellipticCurve(R([1, 2, 1, 2, 6, 4, 1]))
            sage: X.conductor()
            169

        Can't do higher genus yet::

            sage: x = polygen(QQ, 'x')
            sage: C = HyperellipticCurve(x^2+x, x^7+x^2+1)
            sage: C.conductor()
            Traceback (most recent call last):
            ...
            ValueError: curve must be of genus 2

        Sometimes this fails::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(R([0, 1, 0, 0, 0, 1]), R([]));
            sage: C.conductor()
            Traceback (most recent call last):
            ...
            ValueError: even part of the conductor not computed for this curve, use odd_part_only=True to get the odd part

        We can use ``odd_part_only`` in such cases to get the odd part of the conductor at least::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurve(R([0, 1, 0, 0, 0, 1]), R([]));
            sage: C.conductor(odd_part_only=True)
            1

        ::

            sage: R.<x> = PolynomialRing(Rationals())
            sage: C = HyperellipticCurve(R([-30, 0, -37, 0, -15, 0,-2]))
            sage: C.conductor()
            Traceback (most recent call last):
            ...
            ValueError: even part of the conductor not computed for this curve, use odd_part_only=True to get the odd part
            sage: C.conductor(odd_part_only=True)
            15
            sage: C = HyperellipticCurve(R([0, 0, 1, 0, 1]), R([1, 0, 0, 1]))
            sage: C.conductor()
            294
            sage: X = HyperellipticCurve(R([1, 0, 4, 2, 4, 0, 1]))
            sage: X.conductor()
            294

        """
        if not isinstance(self, HyperellipticCurve_g2):
            raise ValueError('curve must be of genus 2')
        rr = self.reduction_type()
        if odd_part_only or (rr.pari_result[1][0][0].sage() != 2 or
                             rr.pari_result[1][1][0].sage() != -1):
            return rr.conductor
        raise ValueError("even part of the conductor not computed for this curve, use odd_part_only=True to get the odd part")

    def lseries(self, prec=53):
        """
        Return the L-series of this hyperelliptic curve of genus 2.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: C = HyperellipticCurve(x^2+x, x^3+x^2+1)
            sage: C.lseries()
            PARI L-function associated to Hyperelliptic Curve
            over Rational Field defined by y^2 + (x^3 + x^2 + 1)*y = x^2 + x
        """
        from sage.lfunctions.pari import LFunction, lfun_genus2
        L = LFunction(lfun_genus2(self), prec=prec)
        L.rename('PARI L-function associated to %s' % self)
        return L
