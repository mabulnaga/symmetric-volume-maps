function [N11,N12,N13,N14,N22,N23,N24,N33,N34,N44] = invertSymmetric4x4_autogenerated(m11,m12,m13,m14,m22,m23,m24,m33,m34,m44)
%INVERTSYMMETRIC4X4_AUTOGENERATED
%    [N11,N12,N13,N14,N22,N23,N24,N33,N34,N44] = INVERTSYMMETRIC4X4_AUTOGENERATED(M11,M12,M13,M14,M22,M23,M24,M33,M34,M44)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    28-Apr-2021 23:34:29

t2 = m12.^2;
t3 = m13.^2;
t4 = m14.^2;
t5 = m23.^2;
t6 = m24.^2;
t7 = m34.^2;
t8 = m11.*m22.*m33.*m44;
t15 = m13.*m14.*m23.*m24.*2.0;
t16 = m12.*m13.*m24.*m34.*2.0;
t17 = m12.*m14.*m23.*m34.*2.0;
t18 = m12.*m14.*m24.*m33.*2.0;
t19 = m13.*m14.*m22.*m34.*2.0;
t20 = m11.*m23.*m24.*m34.*2.0;
t21 = m12.*m13.*m23.*m44.*2.0;
t9 = m11.*m22.*t7;
t10 = m11.*m33.*t6;
t11 = m22.*m33.*t4;
t12 = m11.*m44.*t5;
t13 = m22.*m44.*t3;
t14 = m33.*m44.*t2;
t22 = -t18;
t23 = -t19;
t24 = -t20;
t25 = -t21;
t26 = -t8;
t27 = t3.*t6;
t28 = t4.*t5;
t29 = t2.*t7;
t30 = -t27;
t31 = -t28;
t32 = -t29;
t33 = t9+t10+t11+t12+t13+t14+t15+t16+t17+t22+t23+t24+t25+t26+t30+t31+t32;
t34 = 1.0./t33;
N11 = t34.*(m22.*t7+m33.*t6+m44.*t5-m23.*m24.*m34.*2.0-m22.*m33.*m44);
if nargout > 1
    N12 = -t34.*(m12.*t7-m13.*m24.*m34-m14.*m23.*m34+m14.*m24.*m33+m13.*m23.*m44-m12.*m33.*m44);
end
if nargout > 2
    N13 = -t34.*(m13.*t6-m14.*m23.*m24-m12.*m24.*m34+m14.*m22.*m34+m12.*m23.*m44-m13.*m22.*m44);
end
if nargout > 3
    N14 = -t34.*(m14.*t5-m13.*m23.*m24-m12.*m23.*m34+m12.*m24.*m33+m13.*m22.*m34-m14.*m22.*m33);
end
if nargout > 4
    N22 = t34.*(m11.*t7+m33.*t4+m44.*t3-m13.*m14.*m34.*2.0-m11.*m33.*m44);
end
if nargout > 5
    N23 = -t34.*(m23.*t4-m13.*m14.*m24-m12.*m14.*m34+m11.*m24.*m34+m12.*m13.*m44-m11.*m23.*m44);
end
if nargout > 6
    N24 = -t34.*(m24.*t3-m13.*m14.*m23-m12.*m13.*m34+m12.*m14.*m33+m11.*m23.*m34-m11.*m24.*m33);
end
if nargout > 7
    N33 = t34.*(m11.*t6+m22.*t4+m44.*t2-m12.*m14.*m24.*2.0-m11.*m22.*m44);
end
if nargout > 8
    N34 = -t34.*(m34.*t2-m12.*m13.*m24-m12.*m14.*m23+m13.*m14.*m22+m11.*m23.*m24-m11.*m22.*m34);
end
if nargout > 9
    N44 = t34.*(m11.*t5+m22.*t3+m33.*t2-m12.*m13.*m23.*2.0-m11.*m22.*m33);
end
