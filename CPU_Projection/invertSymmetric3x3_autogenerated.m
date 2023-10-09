function [N11,N12,N13,N22,N23,N33] = invertSymmetric3x3_autogenerated(m11,m12,m13,m22,m23,m33)
%INVERTSYMMETRIC3X3_AUTOGENERATED
%    [N11,N12,N13,N22,N23,N33] = INVERTSYMMETRIC3X3_AUTOGENERATED(M11,M12,M13,M22,M23,M33)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    28-Apr-2021 23:34:29

t2 = m12.^2;
t3 = m13.^2;
t4 = m23.^2;
t5 = m11.*m22.*m33;
t9 = m12.*m13.*m23.*2.0;
t6 = m11.*t4;
t7 = m22.*t3;
t8 = m33.*t2;
t10 = -t9;
t11 = -t5;
t12 = t6+t7+t8+t10+t11;
t13 = 1.0./t12;
N11 = t13.*(t4-m22.*m33);
if nargout > 1
    N12 = -t13.*(m13.*m23-m12.*m33);
end
if nargout > 2
    N13 = -t13.*(m12.*m23-m13.*m22);
end
if nargout > 3
    N22 = t13.*(t3-m11.*m33);
end
if nargout > 4
    N23 = -t13.*(m12.*m13-m11.*m23);
end
if nargout > 5
    N33 = t13.*(t2-m11.*m22);
end