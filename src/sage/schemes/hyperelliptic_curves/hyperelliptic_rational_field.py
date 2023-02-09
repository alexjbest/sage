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

from . import hyperelliptic_generic


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
