function [A, B, C] = ASHRAE(month)
% Returns A [W/m2], B [-], C [-], where:
% Direct normal irradiance In = A * exp(-B / cos(solar zenith))
% Diffuse radition Id = C * In

switch month
    case 1
        A = 1230;
        B = 0.142;
        C = 0.058;
    case 2
        A = 1215;
        B = 0.144;
        C = 0.060;
    case 3
        A = 1186;
        B = 0.156;
        C = 0.097;
    case 4
        A = 1136;
        B = 0.180;
        C = 0.097;
    case 5
        A = 1104;
        B = 0.196;
        C = 0.121;
    case 6
        A = 1088;
        B = 0.205;
        C = 0.134;
    case 7
        A = 1085;
        B = 0.207;
        C = 0.136;
    case 8
        A = 1107;
        B = 0.201;
        C = 0.122;
    case 9
        A = 1152;
        B = 0.177;
        C = 0.092;
    case 10
        A = 1193;
        B = 0.160;
        C = 0.073;
    case 11
        A = 1221;
        B = 0.149;
        C = 0.063;
    case 12
        A = 1234;
        B = 0.142;
        C = 0.057;
end
end