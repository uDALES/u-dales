%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solar Position Algorithm (SPA) for Solar Radiation Application
%                                         
% Measurement & Instrumentation Team Solar Radiation Research Laboratory
% National Renewable Energy Laboratory 1617 Cole Blvd, Golden, CO 80401
%
% Last modified:   2017/03/19   M. Mahooti
%                  2023/10/06   S. Owens  : changed name and added arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spa = solarPosition(time_of_day, longitude, latitude, timezone, elevation)
% clc
% clear
format long g

spa_const

% declare the SPA structure
% enter required input values into SPA structure
spa.year          = time_of_day.Year;%2003;
spa.month         = time_of_day.Month;%10;
spa.day           = time_of_day.Day;%17;
spa.hour          = time_of_day.Hour;%12;
spa.minute        = time_of_day.Minute;%30;
spa.second        = time_of_day.Second;%30;
spa.timezone      = timezone;%-7;
spa.delta_ut1     = 0;
spa.delta_t       = 0;%67;
spa.longitude     = longitude;%-105.1786;
spa.latitude      = latitude;%39.742476;
spa.elevation     = elevation;%1830.14;
spa.pressure      = 1010;%820;
spa.temperature   = 10;%11;
spa.slope         = 0;%30;
spa.azm_rotation  = 0;%-10;
spa.atmos_refract = 0.5667;
spa.function      = SPA_ALL;

%call the SPA calculate function and pass the SPA structure
[result, spa] = spa_calculate(spa);

if (result == 0)  %check for SPA errors
%     %display the results inside the SPA structure
%     fprintf('Julian Day:    %.6f\n',spa.jd);
%     fprintf('L:             %.6e degrees\n',spa.l);
%     fprintf('B:             %.6e degrees\n',spa.b);
%     fprintf('R:             %.6f AU\n',spa.r);
%     fprintf('H:             %.6f degrees\n',spa.h);
%     fprintf('Delta Psi:     %.6e degrees\n',spa.del_psi);
%     fprintf('Delta Epsilon: %.6e degrees\n',spa.del_epsilon);
%     fprintf('Epsilon:       %.6f degrees\n',spa.epsilon);
%     fprintf('Zenith:        %.6f degrees\n',spa.zenith);
%     fprintf('Azimuth:       %.6f degrees\n',spa.azimuth);
%     fprintf('Incidence:     %.6f degrees\n',spa.incidence);
% 
%     min = 60*(spa.sunrise - floor(spa.sunrise));
%     sec = 60*(min - floor(min));
%     fprintf('Sunrise:       %2.2d:%2.2d:%2.2d Local Time\n', ...
%            floor(spa.sunrise), floor(min), floor(sec));
% 
%     min = 60*(spa.sunset - floor(spa.sunset));
%     sec = 60*(min - floor(min));
%     fprintf('Sunset:        %2.2d:%2.2d:%2.2d Local Time\n', ...
%             floor(spa.sunset), floor(min), floor(sec));
else
    fprintf('SPA Error Code: %d\n', result);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output of this program should be:
%
%Julian Day:    2452930.312847
%L:             2.401826e+01 degrees
%B:             -1.011219e-04 degrees
%R:             0.996542 AU
%H:             11.105902 degrees
%Delta Psi:     -3.998404e-03 degrees
%Delta Epsilon: 1.666568e-03 degrees
%Epsilon:       23.440465 degrees
%Zenith:        50.111622 degrees
%Azimuth:       194.340241 degrees
%Incidence:     25.187000 degrees
%Sunrise:       06:12:43 Local Time
%Sunset:        17:20:19 Local Time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limited = limit_degrees(degrees)

degrees = degrees/360;
limited = 360*(degrees-floor(degrees));

if (limited < 0)
    limited = limited + 360;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limited = limit_degrees180pm(degrees)

degrees = degrees/360;
limited = 360*(degrees-floor(degrees));

if(limited < -180)
    limited = limited + 360;
elseif(limited > 180)
    limited = limited - 360;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limited = limit_degrees180(degrees)

degrees = degrees/180;
limited = 180*(degrees-floor(degrees));
if (limited < 0)
    limited = limited + 180;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limited = limit_zero2one(value)

limited = value - floor(value);
if (limited < 0)
    limited = limited + 1;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limited = limit_minutes(minutes)

limited = minutes;

if(limited < -20)
    limited = limited + 1440;
elseif(limited > 20)
    limited = limited - 1440;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = dayfrac_to_local_hr(dayfrac, timezone)

out = 24*limit_zero2one(dayfrac + timezone/24);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = third_order_polynomial(a, b, c, d, x)

out = ((a*x + b)*x + c)*x + d;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = validate_inputs(spa)

spa_const

if ((spa.year        < -2000) || (spa.year        > 6000))
    out = 1;
    return
end
if ((spa.month       < 1    ) || (spa.month       > 12  ))
    out = 2;
    return
end
if ((spa.day         < 1    ) || (spa.day         > 31  ))
    out = 3;
    return
end
if ((spa.hour        < 0    ) || (spa.hour        > 24  ))
    out = 4;
    return
end
if ((spa.minute      < 0    ) || (spa.minute      > 59  ))
    out = 5;
    return
end
if ((spa.second      < 0    ) || (spa.second      >=60  ))
    out = 6;
    return
end
if ((spa.pressure    < 0    ) || (spa.pressure    > 5000))
    out = 12;
    return
end
if ((spa.temperature <= -273) || (spa.temperature > 6000))
    out = 13;
    return
end
if ((spa.delta_ut1   <= -1  ) || (spa.delta_ut1   >= 1  ))
    out = 17;
    return
end
if ((spa.hour        == 24  ) && (spa.minute      > 0   ))
    out = 5;
    return
end
if ((spa.hour        == 24  ) && (spa.second      > 0   ))
    out = 6;
    return
end

if (abs(spa.delta_t)       > 8000    )
    out = 7;
    return
end
if (abs(spa.timezone)      > 18      )
    out = 8;
    return
end
if (abs(spa.longitude)     > 180     )
    out = 9;
    return
end
if (abs(spa.latitude)      > 90      )
    out = 10;
    return
end
if (abs(spa.atmos_refract) > 5       )
    out = 16;
    return
end
if (    spa.elevation      < -6500000)
    out = 11;
    return
end

if ((spa.function == SPA_ZA_INC) || (spa.function == SPA_ALL))
    if (abs(spa.slope)         > 360)
        out = 14;
        return
    end
    if (abs(spa.azm_rotation)  > 360)
        out = 15;
        return
    end
end

out = 0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function julian_day = julian_day(year,month,day,hour,minute,second,dut1,tz)

day_decimal = day+(hour-tz+(minute+(second+dut1)/60)/60)/24;

if(month < 3)
    month = month+12;
    year = year-1;
end

julian_day = floor(365.25*(year+4716))+floor(30.6001*(month+1))+...
             day_decimal-1524.5;

if (julian_day > 2299160)
    a = floor(year/100);
    julian_day = julian_day+(2-a+floor(a/4));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = julian_century(jd)

out = (jd-2451545)/36525;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = julian_ephemeris_day(jd, delta_t)

out = jd+delta_t/86400;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = julian_ephemeris_century(jde)

out = (jde-2451545)/36525;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = julian_ephemeris_millennium(jce)

out = (jce/10);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = earth_periodic_term_summation(terms,count,jme)

spa_const

sum = 0;

for i = 1:count
    sum = sum + (terms(i,TERM_A)*cos(terms(i,TERM_B)+terms(i,TERM_C)*jme));
end

out = sum;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = earth_values(term_sum,count,jme)

sum = 0;

for i = 1:count
    sum = sum + term_sum(i)*(jme^(i-1));
end

sum = sum/1e8;

out = sum;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = earth_heliocentric_longitude(jme)

spa_const

sum = zeros(L_COUNT,1);

for i = 1:L_COUNT
    L_TERM = cell2mat(L_TERMS(i));    
    sum(i) = earth_periodic_term_summation(L_TERM, l_subcount(i), jme);
end

out = limit_degrees(rad2deg(earth_values(sum, L_COUNT, jme)));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = earth_heliocentric_latitude(jme)

spa_const

sum = zeros(B_COUNT,1);

for i = 1:B_COUNT
    B_TERM = cell2mat(B_TERMS(i));
    sum(i) = earth_periodic_term_summation(B_TERM, b_subcount(i), jme);
end

out = rad2deg(earth_values(sum, B_COUNT, jme));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = earth_radius_vector(jme)

spa_const

sum = zeros(R_COUNT,1);

for i = 1:R_COUNT
    R_TERM = cell2mat(R_TERMS(i));
    sum(i) = earth_periodic_term_summation(R_TERM, r_subcount(i), jme);
end

out = earth_values(sum, R_COUNT, jme);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = geocentric_longitude(l)

theta = l+180;

if (theta >= 360)
    theta = theta-360;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = geocentric_latitude(b)

out = -b;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mean_elongation_moon_sun(jce)

out = third_order_polynomial(1/189474,-0.0019142,445267.11148,297.85036,jce);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mean_anomaly_sun(jce)

out = third_order_polynomial(-1/300000,-0.0001603,35999.05034,357.52772,jce);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mean_anomaly_moon(jce)

out = third_order_polynomial(1/56250,0.0086972,477198.867398,134.96298,jce);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = argument_latitude_moon(jce)

out = third_order_polynomial(1/327270,-0.0036825,483202.017538,93.27191,jce);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = ascending_longitude_moon(jce)

out = third_order_polynomial(1/450000,0.0020708,-1934.136261,125.04452,jce);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sum = xy_term_summation(i, x)

spa_const

sum=0;

for j = 1:TERM_Y_COUNT
    sum = sum + x(j)*Y_TERMS(i,j);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [del_psi, del_epsilon] = nutation_longitude_and_obliquity(jce, x)

spa_const

sum_psi = 0;
sum_epsilon = 0;

for i = 1:Y_COUNT
    xy_term_sum = deg2rad(xy_term_summation(i, x));
    sum_psi     = sum_psi + ( (PE_TERMS(i,TERM_PSI_A) + ...
                           jce*PE_TERMS(i,TERM_PSI_B))*sin(xy_term_sum) );
    sum_epsilon = sum_epsilon + ( (PE_TERMS(i,TERM_EPS_C) + ...
                           jce*PE_TERMS(i,TERM_EPS_D))*cos(xy_term_sum) );
end

del_psi     = sum_psi    /36000000;
del_epsilon = sum_epsilon/36000000;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = ecliptic_mean_obliquity(jme)

u = jme/10;

out = 84381.448+u*(-4680.93+u*(-1.55+u*(1999.25+u*(-51.38 + ...
      u*(-249.67+u*(-39.05+u*(7.12+u*(27.87+u*(5.79+u*2.45)))))))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = ecliptic_true_obliquity(delta_epsilon, epsilon0)

out = delta_epsilon + epsilon0/3600;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = aberration_correction(r)

out = -20.4898/(3600*r);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = apparent_sun_longitude(theta, delta_psi, delta_tau)

out = theta + delta_psi + delta_tau;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = greenwich_mean_sidereal_time (jd, jc)

out = limit_degrees(280.46061837+360.98564736629*(jd-2451545)+ ...
      jc*jc*(0.000387933 - jc/38710000));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = greenwich_sidereal_time (nu0, delta_psi, epsilon)

out = nu0 + delta_psi*cos(deg2rad(epsilon));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = geocentric_right_ascension(lamda, epsilon, beta)

lamda_rad   = deg2rad(lamda);
epsilon_rad = deg2rad(epsilon);

out = limit_degrees(rad2deg(atan2(sin(lamda_rad)*cos(epsilon_rad) - ...
      tan(deg2rad(beta))*sin(epsilon_rad), cos(lamda_rad))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = geocentric_declination(beta, epsilon, lamda)

beta_rad    = deg2rad(beta);
epsilon_rad = deg2rad(epsilon);

out = rad2deg(asin(sin(beta_rad)*cos(epsilon_rad) + ...
      cos(beta_rad)*sin(epsilon_rad)*sin(deg2rad(lamda))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = observer_hour_angle(nu, longitude, alpha_deg)

out = limit_degrees(nu + longitude - alpha_deg);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sun_equatorial_horizontal_parallax(r)

out = 8.794/(3600*r);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta_alpha, delta_prime] = ...
right_ascension_parallax_and_topocentric_dec(latitude,elevation,xi,h,delta)

lat_rad   = deg2rad(latitude);
xi_rad    = deg2rad(xi);
h_rad     = deg2rad(h);
delta_rad = deg2rad(delta);
u = atan(0.99664719 * tan(lat_rad));
y = 0.99664719 * sin(u) + elevation*sin(lat_rad)/6378140;
x =              cos(u) + elevation*cos(lat_rad)/6378140;

delta_alpha_rad = atan2(-x*sin(xi_rad)*sin(h_rad), ...
                        cos(delta_rad)-x*sin(xi_rad)*cos(h_rad));

delta_prime = rad2deg(atan2((sin(delta_rad)-y*sin(xi_rad))* ...
            cos(delta_alpha_rad),cos(delta_rad)-x*sin(xi_rad)*cos(h_rad)));

delta_alpha = rad2deg(delta_alpha_rad);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_right_ascension(alpha_deg, delta_alpha)

out = alpha_deg + delta_alpha;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_local_hour_angle(h, delta_alpha)

out = h - delta_alpha;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_elevation_angle(latitude, delta_prime, h_prime)

lat_rad         = deg2rad(latitude);
delta_prime_rad = deg2rad(delta_prime);

out = rad2deg(asin(sin(lat_rad)*sin(delta_prime_rad) + ...
              cos(lat_rad)*cos(delta_prime_rad) * cos(deg2rad(h_prime))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function del_e = atmospheric_refraction_correction(pressure,temperature,...
	                                                     atmos_refract, e0)
spa_const

del_e = 0;

if(e0 >= -1*(SUN_RADIUS + atmos_refract))
    del_e = (pressure/1010)*(283/(273+temperature))* ...
             1.02/(60*tan(deg2rad(e0+10.3/(e0+5.11))));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_elevation_angle_corrected(e0, delta_e)

out = e0 + delta_e;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_zenith_angle(e)

out = 90.0 - e;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_azimuth_angle_astro(h_prime, latitude, delta_prime)

h_prime_rad = deg2rad(h_prime);
lat_rad     = deg2rad(latitude);

out = limit_degrees(rad2deg(atan2(sin(h_prime_rad), ...
 cos(h_prime_rad)*sin(lat_rad) - tan(deg2rad(delta_prime))*cos(lat_rad))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = topocentric_azimuth_angle(azimuth_astro)

out = limit_degrees(azimuth_astro + 180);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = surface_incidence_angle(zenith,azimuth_astro,azm_rotation,slope)

zenith_rad = deg2rad(zenith);
slope_rad  = deg2rad(slope);

out = rad2deg(acos(cos(zenith_rad)*cos(slope_rad) + ...
 sin(slope_rad)*sin(zenith_rad)*cos(deg2rad(azimuth_astro-azm_rotation))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sun_mean_longitude(jme)

out = limit_degrees(280.4664567+jme*(360007.6982779+jme*(0.03032028+ ...
                    jme*(1/49931+jme*(-1/15300+jme*(-1/2000000))))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eot(m, alpha, del_psi, epsilon)

out = limit_minutes(4*(m-0.0057183-alpha+del_psi*cos(deg2rad(epsilon))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = approx_sun_transit_time(alpha_zero, longitude, nu)

out = (alpha_zero-longitude-nu)/360;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h0 = sun_hour_angle_at_rise_set(latitude, delta_zero, h0_prime)

h0             = -99999;
latitude_rad   = deg2rad(latitude);
delta_zero_rad = deg2rad(delta_zero);
argument       = (sin(deg2rad(h0_prime))-sin(latitude_rad)* ...
              sin(delta_zero_rad))/(cos(latitude_rad)*cos(delta_zero_rad));
if (abs(argument) <= 1)
    h0 = limit_degrees180(rad2deg(acos(argument)));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m_rts = approx_sun_rise_and_set(m_rts, h0)

spa_const

h0_dfrac = h0/360;

m_rts(SUN_RISE)    = limit_zero2one(m_rts(SUN_TRANSIT) - h0_dfrac);
m_rts(SUN_SET)     = limit_zero2one(m_rts(SUN_TRANSIT) + h0_dfrac);
m_rts(SUN_TRANSIT) = limit_zero2one(m_rts(SUN_TRANSIT));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rts_alpha_delta_prime(ad, n)

spa_const

a = ad(JD_ZERO) - ad(JD_MINUS);
b = ad(JD_PLUS) - ad(JD_ZERO);

if (abs(a) >= 2)
    a = limit_zero2one(a);
end
if (abs(b) >= 2)
    b = limit_zero2one(b);
end

out = ad(JD_ZERO) + n * (a + b + (b-a)*n)/2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rts_sun_altitude(latitude, delta_prime, h_prime)

latitude_rad    = deg2rad(latitude);
delta_prime_rad = deg2rad(delta_prime);

out = rad2deg(asin(sin(latitude_rad)*sin(delta_prime_rad) + ...
            cos(latitude_rad)*cos(delta_prime_rad)*cos(deg2rad(h_prime))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sun_rise_and_set(m_rts, h_rts, delta_prime, latitude, ...
                                h_prime, h0_prime, sun)

out = m_rts(sun)+(h_rts(sun)-h0_prime)/(360*cos(deg2rad(delta_prime(sun)))...
                *cos(deg2rad(latitude))*sin(deg2rad(h_prime(sun))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate required SPA parameters to get the right ascension (alpha) and
% declination (delta) Note: JD must be already calculated and in structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spa = calculate_geocentric_sun_right_ascension_and_declination(spa)

spa_const

x = zeros(TERM_X_COUNT,1);

spa.jc = julian_century(spa.jd);

spa.jde = julian_ephemeris_day(spa.jd, spa.delta_t);
spa.jce = julian_ephemeris_century(spa.jde);
spa.jme = julian_ephemeris_millennium(spa.jce);

spa.l = earth_heliocentric_longitude(spa.jme);
spa.b = earth_heliocentric_latitude(spa.jme);
spa.r = earth_radius_vector(spa.jme);

spa.theta = geocentric_longitude(spa.l);
spa.beta  = geocentric_latitude(spa.b);

spa.x0 = mean_elongation_moon_sun(spa.jce);
x(TERM_X0) = spa.x0;
spa.x1 = mean_anomaly_sun(spa.jce);
x(TERM_X1) = spa.x1;
spa.x2 = mean_anomaly_moon(spa.jce);
x(TERM_X2) = spa.x2;
spa.x3 = argument_latitude_moon(spa.jce);
x(TERM_X3) = spa.x3;
spa.x4 = ascending_longitude_moon(spa.jce);
x(TERM_X4) = spa.x4;

[spa.del_psi,spa.del_epsilon] = nutation_longitude_and_obliquity(spa.jce,x);

spa.epsilon0 = ecliptic_mean_obliquity(spa.jme);
spa.epsilon  = ecliptic_true_obliquity(spa.del_epsilon, spa.epsilon0);

spa.del_tau = aberration_correction(spa.r);
spa.lamda   = apparent_sun_longitude(spa.theta, spa.del_psi,spa.del_tau);
spa.nu0     = greenwich_mean_sidereal_time (spa.jd, spa.jc);
spa.nu      = greenwich_sidereal_time (spa.nu0, spa.del_psi,spa.epsilon);

spa.alpha = geocentric_right_ascension(spa.lamda, spa.epsilon,spa.beta);
spa.delta = geocentric_declination(spa.beta, spa.epsilon, spa.lamda);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spa = calculate_eot_and_sun_rise_transit_set(spa)

spa_const

alpha = zeros(JD_COUNT,1);
delta = zeros(JD_COUNT,1);
m_rts = zeros(SUN_COUNT,1);
nu_rts = zeros(SUN_COUNT,1);
h_rts = zeros(SUN_COUNT,1);
alpha_prime = zeros(SUN_COUNT,1);
delta_prime = zeros(SUN_COUNT,1);
h_prime = zeros(SUN_COUNT,1);
h0_prime = -1*(SUN_RADIUS + spa.atmos_refract);

sun_rts = spa;
m       = sun_mean_longitude(spa.jme);
spa.eot = eot(m, spa.alpha, spa.del_psi, spa.epsilon);
sun_rts.hour = 0;
sun_rts.minute = 0;
sun_rts.second = 0;
sun_rts.delta_ut1 = 0;
sun_rts.timezone = 0;

sun_rts.jd = julian_day(sun_rts.year,sun_rts.month,sun_rts.day,sun_rts.hour, ...
         sun_rts.minute,sun_rts.second,sun_rts.delta_ut1,sun_rts.timezone);

sun_rts = calculate_geocentric_sun_right_ascension_and_declination(sun_rts);
nu = sun_rts.nu;

sun_rts.delta_t = 0;
sun_rts.jd = sun_rts.jd-1;
for i = 1:JD_COUNT
    sun_rts = calculate_geocentric_sun_right_ascension_and_declination(sun_rts);
    alpha(i) = sun_rts.alpha;
    delta(i) = sun_rts.delta;
    sun_rts.jd = sun_rts.jd+1;
end

m_rts(SUN_TRANSIT) = approx_sun_transit_time(alpha(JD_ZERO),spa.longitude,nu);
h0 = sun_hour_angle_at_rise_set(spa.latitude,delta(JD_ZERO),h0_prime);

if (h0 >= 0)
    m_rts = approx_sun_rise_and_set(m_rts, h0);
    for i = 1:SUN_COUNT
        nu_rts(i)      = nu + 360.985647*m_rts(i);
        n              = m_rts(i) + spa.delta_t/86400;
        alpha_prime(i) = rts_alpha_delta_prime(alpha,n);
        delta_prime(i) = rts_alpha_delta_prime(delta,n);
        h_prime(i)     = limit_degrees180pm(nu_rts(i)+spa.longitude-alpha_prime(i));
        h_rts(i)       = rts_sun_altitude(spa.latitude,delta_prime(i),h_prime(i));
    end
    
    spa.srha = h_prime(SUN_RISE);
    spa.ssha = h_prime(SUN_SET);
    spa.sta  = h_rts(SUN_TRANSIT);
    
    spa.suntransit = dayfrac_to_local_hr(m_rts(SUN_TRANSIT)- ...
                                    h_prime(SUN_TRANSIT)/360,spa.timezone);
    
    spa.sunrise = dayfrac_to_local_hr(sun_rise_and_set(m_rts,h_rts, ...
         delta_prime,spa.latitude,h_prime,h0_prime,SUN_RISE),spa.timezone);
    
    spa.sunset = dayfrac_to_local_hr(sun_rise_and_set(m_rts,h_rts, ...
          delta_prime,spa.latitude,h_prime,h0_prime,SUN_SET),spa.timezone);
    
else
    spa.srha = -99999;
    spa.ssha = -99999;
    spa.sta = -99999;
    spa.suntransit = -99999;
    spa.sunrise = -99999;
    spa.sunset = -99999;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate all SPA parameters and put into structure
% Note: All inputs values (listed in header file) must already be in structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, spa] = spa_calculate(spa)

spa_const

result = validate_inputs(spa);

if (result == 0)
    spa.jd = julian_day(spa.year,spa.month,spa.day,spa.hour,spa.minute, ...
                                  spa.second, spa.delta_ut1, spa.timezone);
    
    spa = calculate_geocentric_sun_right_ascension_and_declination(spa);
    
    spa.h  = observer_hour_angle(spa.nu, spa.longitude, spa.alpha);
    spa.xi = sun_equatorial_horizontal_parallax(spa.r);

    [spa.del_alpha, spa.delta_prime] = ...
    right_ascension_parallax_and_topocentric_dec(spa.latitude, ...
                                     spa.elevation,spa.xi,spa.h,spa.delta);

    spa.alpha_prime = topocentric_right_ascension(spa.alpha, spa.del_alpha);
    spa.h_prime     = topocentric_local_hour_angle(spa.h, spa.del_alpha);

    spa.e0    = topocentric_elevation_angle(spa.latitude, spa.delta_prime, spa.h_prime);
    spa.del_e = atmospheric_refraction_correction(spa.pressure, ...
                                spa.temperature,spa.atmos_refract, spa.e0);
    spa.e     = topocentric_elevation_angle_corrected(spa.e0, spa.del_e);

    spa.zenith        = topocentric_zenith_angle(spa.e);
    spa.azimuth_astro = topocentric_azimuth_angle_astro(spa.h_prime, ...
                                             spa.latitude,spa.delta_prime);
    spa.azimuth       = topocentric_azimuth_angle(spa.azimuth_astro);

    if ((spa.function == SPA_ZA_INC) || (spa.function == SPA_ALL))
        spa.incidence  = surface_incidence_angle(spa.zenith, ...
                             spa.azimuth_astro,spa.azm_rotation,spa.slope);
    end
    if ((spa.function == SPA_ZA_RTS) || (spa.function == SPA_ALL))
        spa = calculate_eot_and_sun_rise_transit_set(spa);
    end
end

end

