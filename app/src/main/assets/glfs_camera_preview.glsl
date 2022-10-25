
#extension GL_OES_EGL_image_external: require

precision mediump float;

uniform samplerExternalOES sTexture;

varying vec2 vTextureCoord;

//
// Helpers
//
const float PI = 3.141592653589793;
const float FLT_MAX = 3.402823466e+38;

float degreesToRadians(float x) {
    return x * PI / 180.0;
}

float radiansToDegrees(float x) {
    return x * 180.0 / PI;
}

float atan2(float y, float x) {
    float s = (abs(x) > abs(y)) ? 1.0 : 0.0;
    return mix(PI / 2.0 - atan(x, y), atan(y, x), s);
}

float cbrt(float x) {
    return sign(x) * pow(abs(x), 1.0 / 3.0);
}

float square(float x) {
    return x * x;
}

float cube(float x) {
    return x * x * x;
}

float log10(in float n) {
    const float kLogBase10 = 1.0 / log2(10.0);
    return log2(n) * kLogBase10;
}

vec3 crosstalk(vec3 x, float a) {
    float b = 1.0 - 2.0 * a;
    mat3 M = mat3(b, a, a, a, b, a, a, a, b);
    return M * x;
}

vec3 inverse_crosstalk(vec3 x, float a) {
    float b = 1.0 - a;
    float c = 1.0 - 3.0 * a;
    mat3 M = mat3(b, -a, -a, -a, b, -a, -a, -a, b) / c;
    return M * x;
}

//
// Gamut conversion
//
const mat3 M_SRGB_TO_XYZ = mat3(
0.4123907992659595, 0.2126390058715103,  0.01933081871559185,
0.3575843393838779, 0.7151686787677559,  0.1191947797946259,
0.1804807884018343, 0.07219231536073372, 0.9505321522496606
);
const mat3 M_XYZ_TO_SRGB = mat3(
3.240969941904521, -0.9692436362808798,  0.05563007969699361,
-1.537383177570093,  1.87596750150772,   -0.2039769588889765,
-0.4986107602930033, 0.04155505740717561, 1.056971514242878
);

const mat3 M_DISPLAY_P3_TO_XYZ = mat3(
0.4865709486482162, 0.2289745640697488, 0.0,
0.2656676931690929, 0.6917385218365062, 0.04511338185890257,
0.1982172852343625, 0.079286914093745,  1.04394436890097500
);
const mat3 M_XYZ_TO_DISPLAY_P3 = mat3(
2.493496911941424, -0.829488969561575,   0.03584583024378433,
-0.9313836179191236, 1.762664060318346,  -0.07617238926804171,
-0.4027107844507168, 0.02362468584194359, 0.9568845240076873
);

const mat3 M_BT2020_TO_XYZ = mat3(
0.6369580483012913, 0.262700212011267,   0.0,
0.1446169035862083, 0.677998071518871,   0.0280726930490875,
0.168880975164172,  0.05930171646986194, 1.06098505771079
);
const mat3 M_XYZ_TO_BT2020 = mat3(
1.716651187971267, -0.666684351832489,   0.01763985744531091,
-0.3556707837763924, 1.616481236634939,  -0.04277061325780865,
-0.2533662813736598, 0.01576854581391113, 0.942103121235474
);

vec3 xyzToOklab(vec3 c) {
    float l = 0.8189330101 * c.x + 0.3618667424 * c.y - 0.1288597137 * c.z;
    float m = 0.0329845436 * c.x + 0.9293118715 * c.y + 0.0361456387 * c.z;
    float s = 0.0482003018 * c.x + 0.2643662691 * c.y + 0.6338517070 * c.z;

    float l_ = cbrt(l);
    float m_ = cbrt(m);
    float s_ = cbrt(s);

    float L = 0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_;
    float a = 1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_;
    float b = 0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_;

    return vec3(L, a, b);
}

vec3 oklabToXyz(vec3 c) {
    float L = c.x;
    float a = c.y;
    float b = c.z;

    float l_ = L + 0.3963377774 * a + 0.2158037573 * b;
    float m_ = L - 0.1055613458 * a - 0.0638541728 * b;
    float s_ = L - 0.0894841775 * a - 1.2914855480 * b;

    float l = l_ * l_ * l_;
    float m = m_ * m_ * m_;
    float s = s_ * s_ * s_;

    return vec3(
    +1.2270138511 * l - 0.5577999807 * m + 0.2812561490 * s,
    -0.0405801784 * l + 1.1122568696 * m - 0.0716766787 * s,
    -0.0763812845 * l - 0.4214819784 * m + 1.5861632204 * s
    );
}

vec3 linearSrgbToXyz(vec3 c) {
    return M_SRGB_TO_XYZ * c;
}

vec3 xyzToLinearSrgb(vec3 c) {
    return M_XYZ_TO_SRGB * c;
}

vec3 linearDisplayP3ToXyz(vec3 c) {
    return M_DISPLAY_P3_TO_XYZ * c;
}

vec3 xyzToLinearDisplayP3(vec3 c) {
    return M_XYZ_TO_DISPLAY_P3 * c;
}

vec3 linearBT2020ToXyz(vec3 c) {
    return M_BT2020_TO_XYZ * c;
}

vec3 xyzToLinearBT2020(vec3 c) {
    return M_XYZ_TO_BT2020 * c;
}

//vec3 linearSrgbToOklab(vec3 c) {
//    return xyzToOklab(linearSrgbToXyz(c));
//}
//
//vec3 oklabToLinearSrgb(vec3 c) {
//    return xyzToLinearSrgb(oklabToXyz(c));
//}

vec3 linearDisplayP3ToOklab(vec3 c) {
    return xyzToOklab(linearDisplayP3ToXyz(c));
}

vec3 oklabToLinearDisplayP3(vec3 c) {
    return xyzToLinearDisplayP3(oklabToXyz(c));
}

vec3 linearBT2020ToOklab(vec3 c) {
    return xyzToOklab(linearBT2020ToXyz(c));
}

vec3 oklabToLinearBT2020(vec3 c) {
    return xyzToLinearBT2020(oklabToXyz(c));
}

//
// OETF, EOTF, OOTF
//

float srgb_to_linear(float srgb) {
    return srgb <= 0.04045 ? srgb / 12.92 : pow((srgb + 0.055) / 1.055, 2.4);
}

float srgb_from_linear(float linear) {
    return linear <= 0.0031308 ? 12.92 * linear : 1.055 * pow(linear, 1.0 / 2.4) - 0.055;
}

const highp float m1 = (2610.0 / 4096.0) / 4.0;
const highp float m2 = (2523.0 / 4096.0) * 128.0;
const highp float c1 = (3424.0 / 4096.0);
const highp float c2 = (2413.0 / 4096.0) * 32.0;
const highp float c3 = (2392.0 / 4096.0) * 32.0;

highp float pq_to_linear(highp float pq) {
    highp float pq_pow_inv_m2 = pow(pq, 1.0 / m2);
    return pow(max(0.0, pq_pow_inv_m2 - c1) / (c2 - c3 * pq_pow_inv_m2), 1.0 / m1);
}

highp float pq_from_linear(highp float linear) {
    highp float linear_pow_m1 = pow(linear, m1);
    return pow((c1 + c2 * linear_pow_m1) / (1.0 + c3 * linear_pow_m1), m2);
}

// Note: it is perhaps debatable whether “linear” for HLG should be scene light
// or display light. Here, it is implemented in terms of scene light.
const float kHlgA = 0.17883277;
const float kHlgB = 0.28466892;
const float kHlgC = 0.55991073;

float hlg_to_linear(float hlg) {
    return hlg <= 0.5 ? hlg * hlg / 3.0 : (exp((hlg - kHlgC) / kHlgA) + kHlgB) / 12.0;
}
float hlg_from_linear(float linear) {
    return linear <= 1.0 / 12.0 ? sqrt(3.0 * linear) : kHlgA * log(12.0 * linear - kHlgB) + kHlgC;
}

// HLG EOTF
// Generate EOTF that converts signal values to relative scene light, both normalized to [0, 1]
highp float EOTF_channel(const highp float channel) {
    const highp float a = 0.17883277;
    const highp float b = 0.28466892; // 1.-4.*a;
    const highp float c = 0.55991073; // 0.5-a*log(4.*a);
    return channel <= 0.5 ? channel * channel / 3.0 : (exp((channel - c) / a) + b) / 12.0;
}

vec3 EOTF(const highp vec3 color) {
    return clamp(vec3(EOTF_channel(color.r), EOTF_channel(color.g), EOTF_channel(color.b)), 0.0, 1.0);
}

highp float HLG_inv_OETF(const highp float channel) {
    return EOTF_channel(channel);
}

highp float HLG_OOTF(const highp float channel, highp float lambada) {
    return channel * pow(channel, lambada - 1.0);
}

highp vec3 HLG_OOTF_3f(highp vec3 color, highp float lambada) {
    highp float y = color.r * 0.262700212011267 + color.g * 0.677998071518871 + color.b * 0.05930171646986194;
    return color * pow(y, lambada - 1.0);
}

// sRGB OETF
// Generate OETF that converts relative display light to signal values, both normalized to [0, 1]
float OETF_sRGB(const float linear) {
    return linear <= 0.0031308 ? linear * 12.92 : (pow(linear, 1.0 / 2.4) * 1.055) - 0.055;
}

vec3 OETF_sRGB(const vec3 linear) {
    return clamp(vec3(OETF_sRGB(linear.r), OETF_sRGB(linear.g), OETF_sRGB(linear.b)), 0.0, 1.0);
}

vec3 OETF(const vec3 linear) {
    return sign(linear.rgb) * OETF_sRGB(abs(linear.rgb));
}

const highp float pq_m1 = 0.1593017578125; // (2610.0 / 4096.0) / 4.0;
const highp float pq_m2 = 78.84375; // (2523.0 / 4096.0) * 128.0;
const highp float pq_c1 = 0.8359375; // (3424.0 / 4096.0);
const highp float pq_c2 = 18.8515625; // (2413.0 / 4096.0) * 32.0;
const highp float pq_c3 = 18.6875; // (2392.0 / 4096.0) * 32.0;

const float pq_C = 10000.0;

// Converts from the non-linear perceptually quantized space to linear cd/m^2
// Note that this is in float, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex
// sections of SMPTE ST 2084-2014
highp float ST2084_2_Y(const highp float N) {
    // Note that this does NOT handle any of the signal range
    // considerations from 2084 - this assumes full range (0 - 1)
    highp float Np = pow( N, 1.0 / pq_m2 );
    highp float L = Np - pq_c1;
    if ( L < 0.0 ) {
        L = 0.0;
    }
    L = L / ( pq_c2 - pq_c3 * Np );
    L = pow( L, 1.0 / pq_m1 );
    return L * pq_C; // returns cd/m^2
}

// Converts from linear cd/m^2 to the non-linear perceptually quantized space
// Note that this is in float, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex
// sections of SMPTE ST 2084-2014
highp float Y_2_ST2084( float C )
//pq_r
{
    // Note that this does NOT handle any of the signal range
    // considerations from 2084 - this returns full range (0 - 1)
    highp float L = C / pq_C;
    highp float Lm = pow( L, pq_m1 );
    highp float N = ( pq_c1 + pq_c2 * Lm ) / ( 1.0 + pq_c3 * Lm );
    N = pow( N, pq_m2 );
    return N;
}

vec3 Y_2_ST2084_f3( vec3 color)
{// converts from linear cd/m^2 to PQ code values
    vec3 outColor;
    outColor.x = Y_2_ST2084( color.x );
    outColor.y = Y_2_ST2084( color.y );
    outColor.z = Y_2_ST2084( color.z );

    return outColor;
}

vec3 ST2084_2_Y_f3( vec3 color)
{
// converts from PQ code values to linear cd/m^2
    vec3 outColor;
    outColor.x = ST2084_2_Y( color.x);
    outColor.y = ST2084_2_Y( color.y);
    outColor.z = ST2084_2_Y( color.z);
    return outColor;
}

// Conversion of PQ signal to HLG, as detailed in Section 7 of ITU-R BT.2390-0
vec3 ST2084_2_HLG_1000nits_f3( vec3 PQ)
{
    // ST.2084 EOTF (non-linear PQ to display light)
    vec3 displayLinear = ST2084_2_Y_f3( PQ) * pq_C;

    // HLG Inverse EOTF (i.e. HLG inverse OOTF followed by the HLG OETF)
    // HLG Inverse OOTF (display linear to scene linear)
    float Y_d = 0.2627*displayLinear.r + 0.6780*displayLinear.g + 0.0593*displayLinear.b;
    const float L_w = 1000.;
    const float L_b = 0.;
    const float alpha = (L_w-L_b);
    const float beta = L_b;
    const float gamma = 1.2;

    vec3 sceneLinear;
    if (Y_d == 0.) {
        /* This case is to protect against pow(0,-N)=Inf error. The ITU document
        does not offer a recommendation for this corner case. There may be a
        better way to handle this, but for now, this works.
        */
        sceneLinear.r = 0.;
        sceneLinear.g = 0.;
        sceneLinear.b = 0.;
    } else {
        sceneLinear.r = pow( (Y_d-beta)/alpha, (1.-gamma)/gamma) * ((displayLinear.r-beta)/alpha);
        sceneLinear.g = pow( (Y_d-beta)/alpha, (1.-gamma)/gamma) * ((displayLinear.g-beta)/alpha);
        sceneLinear.b = pow( (Y_d-beta)/alpha, (1.-gamma)/gamma) * ((displayLinear.b-beta)/alpha);
    }

    // HLG OETF (scene linear to non-linear signal value)
    const float a = 0.17883277;
    const float b = 0.28466892; // 1.-4.*a;
    const float c = 0.55991073; // 0.5-a*log(4.*a);

    vec3 HLG;
    if (sceneLinear.r <= 1./12.) {
        HLG.r = sqrt(3.*sceneLinear.r);
    } else {
        HLG.r = a*log(12.*sceneLinear.r-b)+c;
    }
    if (sceneLinear.g <= 1./12.) {
        HLG.g = sqrt(3.*sceneLinear.g);
    } else {
        HLG.g = a*log(12.*sceneLinear.g-b)+c;
    }
    if (sceneLinear.b <= 1./12.) {
        HLG.b = sqrt(3.*sceneLinear.b);
    } else {
        HLG.b = a*log(12.*sceneLinear.b-b)+c;
    }

    return HLG;
}



//
// Gamut mapping
//
// Finds the maximum saturation possible for a given hue that fits in sRGB
// Saturation here is defined as S = C/L
// a and b must be normalized so a^2 + b^2 == 1
float compute_max_saturation(float a, float b) {
    // Max saturation will be when one of r, g or b goes below zero.

    // Select different coefficients depending on which component goes below zero first
    float k0, k1, k2, k3, k4, wl, wm, ws;

    if (-1.88170328 * a - 0.80936493 * b > 1.0) {
        // Red component
        k0 = +1.19086277;
        k1 = +1.76576728;
        k2 = +0.59662641;
        k3 = +0.75515197;
        k4 = +0.56771245;
        wl = +4.0767416621;
        wm = -3.3077115913;
        ws = +0.2309699292;
    } else if (1.81444104 * a - 1.19445276 * b > 1.0) {
        // Green component
        k0 = +0.73956515;
        k1 = -0.45954404;
        k2 = +0.08285427;
        k3 = +0.12541070;
        k4 = +0.14503204;
        wl = -1.2681437731;
        wm = +2.6097574011;
        ws = -0.3413193965;
    } else {
        // Blue component
        k0 = +1.35733652;
        k1 = -0.00915799;
        k2 = -1.15130210;
        k3 = -0.50559606;
        k4 = +0.00692167;
        wl = -0.0041960863;
        wm = -0.7034186147;
        ws = +1.7076147010;
    }

    // Approximate max saturation using a polynomial:
    float S = k0 + k1 * a + k2 * b + k3 * a * a + k4 * a * b;

    // Do one step Halley's method to get closer
    // this gives an error less than 10e6, except for some blue hues where the dS/dh is close to infinite
    // this should be sufficient for most applications, otherwise do two/three steps

    float k_l = +0.3963377774 * a + 0.2158037573 * b;
    float k_m = -0.1055613458 * a - 0.0638541728 * b;
    float k_s = -0.0894841775 * a - 1.2914855480 * b;

    {
        float l_ = 1.0 + S * k_l;
        float m_ = 1.0 + S * k_m;
        float s_ = 1.0 + S * k_s;

        float l = l_ * l_ * l_;
        float m = m_ * m_ * m_;
        float s = s_ * s_ * s_;

        float l_dS = 3.0 * k_l * l_ * l_;
        float m_dS = 3.0 * k_m * m_ * m_;
        float s_dS = 3.0 * k_s * s_ * s_;

        float l_dS2 = 6.0 * k_l * k_l * l_;
        float m_dS2 = 6.0 * k_m * k_m * m_;
        float s_dS2 = 6.0 * k_s * k_s * s_;

        float f = wl * l + wm * m + ws * s;
        float f1 = wl * l_dS + wm * m_dS + ws * s_dS;
        float f2 = wl * l_dS2 + wm * m_dS2 + ws * s_dS2;

        S = S - f * f1 / (f1 * f1 - 0.5 * f * f2);
    }

    return S;
}

// finds L_cusp and C_cusp for a given hue
// a and b must be normalized so a^2 + b^2 == 1
vec2 find_cusp(float a, float b) {
    // First, find the maximum saturation (saturation S = C/L)
    float S_cusp = compute_max_saturation(a, b);

    // Convert to linear sRGB to find the first point where at least one of r,g or b >= 1:
    vec3 rgb_at_max = oklabToLinearDisplayP3(vec3(1.0, S_cusp * a, S_cusp * b));
    float L_cusp = cbrt(1.0 / max(max(rgb_at_max.r, rgb_at_max.g), rgb_at_max.b));
    float C_cusp = L_cusp * S_cusp;

    return vec2(L_cusp, C_cusp);
}

// Finds intersection of the line defined by
// L = L0 * (1 - t) + t * L1;
// C = t * C1;
// a and b must be normalized so a^2 + b^2 == 1
float find_gamut_intersection(float a, float b, float L1, float C1, float L0) {
    // Find the cusp of the gamut triangle
    vec2 cusp = find_cusp(a, b);
    float cuspL = cusp.x;
    float cuspC = cusp.y;

    // Find the intersection for upper and lower half seprately
    float t;
    if (((L1 - L0) * cuspC - (cuspL - L0) * C1) <= 0.0) {
        // Lower half

        t = cuspC * L0 / (C1 * cuspL + cuspC * (L0 - L1));
    } else {
        // Upper half

        // First intersect with triangle
        t = cuspC * (L0 - 1.0) / (C1 * (cuspL - 1.0) + cuspC * (L0 - L1));

        // Then one step Halley's method
        {
            float dL = L1 - L0;
            float dC = C1;

            float k_l = +0.3963377774 * a + 0.2158037573 * b;
            float k_m = -0.1055613458 * a - 0.0638541728 * b;
            float k_s = -0.0894841775 * a - 1.2914855480 * b;

            float l_dt = dL + dC * k_l;
            float m_dt = dL + dC * k_m;
            float s_dt = dL + dC * k_s;

            // If higher accuracy is required, 2 or 3 iterations of the following block can be used:
            {
                float L = L0 * (1.0 - t) + t * L1;
                float C = t * C1;

                float l_ = L + C * k_l;
                float m_ = L + C * k_m;
                float s_ = L + C * k_s;

                float l = l_ * l_ * l_;
                float m = m_ * m_ * m_;
                float s = s_ * s_ * s_;

                float ldt = 3.0 * l_dt * l_ * l_;
                float mdt = 3.0 * m_dt * m_ * m_;
                float sdt = 3.0 * s_dt * s_ * s_;

                float ldt2 = 6.0 * l_dt * l_dt * l_;
                float mdt2 = 6.0 * m_dt * m_dt * m_;
                float sdt2 = 6.0 * s_dt * s_dt * s_;

                float r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s - 1.0;
                float r1 = 4.0767416621 * ldt - 3.3077115913 * mdt + 0.2309699292 * sdt;
                float r2 = 4.0767416621 * ldt2 - 3.3077115913 * mdt2 + 0.2309699292 * sdt2;

                float u_r = r1 / (r1 * r1 - 0.5 * r * r2);
                float t_r = -r * u_r;

                float g = -1.2681437731 * l + 2.6097574011 * m - 0.3413193965 * s - 1.0;
                float g1 = -1.2681437731 * ldt + 2.6097574011 * mdt - 0.3413193965 * sdt;
                float g2 = -1.2681437731 * ldt2 + 2.6097574011 * mdt2 - 0.3413193965 * sdt2;

                float u_g = g1 / (g1 * g1 - 0.5 * g * g2);
                float t_g = -g * u_g;

                float b = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s - 1.0;
                float b1 = -0.0041960863 * ldt - 0.7034186147 * mdt + 1.7076147010 * sdt;
                float b2 = -0.0041960863 * ldt2 - 0.7034186147 * mdt2 + 1.7076147010 * sdt2;

                float u_b = b1 / (b1 * b1 - 0.5 * b * b2);
                float t_b = -b * u_b;

                t_r = u_r >= 0.0 ? t_r : FLT_MAX;
                t_g = u_g >= 0.0 ? t_g : FLT_MAX;
                t_b = u_b >= 0.0 ? t_b : FLT_MAX;

                t += min(t_r, min(t_g, t_b));
            }
        }
    }

    return t;
}

vec3 gamut_clip_preserve_lightness(vec3 rgb) {
    if (rgb.r < 1.0 && rgb.g < 1.0 && rgb.b < 1.0 && rgb.r > 0.0 && rgb.g > 0.0 && rgb.b > 0.0)
    return rgb;

    vec3 lab = linearBT2020ToOklab(rgb);

    float L = lab.x;
    float eps = 0.00001;
    float C = max(eps, sqrt(lab.y * lab.y + lab.z * lab.z));
    float a_ = lab.y / C;
    float b_ = lab.z / C;

    float L0 = clamp(L, 0.0, 1.0);

    float t = find_gamut_intersection(a_, b_, L, C, L0);
    float L_clipped = L0 * (1.0 - t) + t * L;
    float C_clipped = t * C;

    return oklabToLinearDisplayP3(vec3(L_clipped, C_clipped * a_, C_clipped * b_));
}

vec3 gamut_clip_project_to_0_5(vec3 rgb) {
    if (rgb.r < 1.0 && rgb.g < 1.0 && rgb.b < 1.0 && rgb.r > 0.0 && rgb.g > 0.0 && rgb.b > 0.0)
    return rgb;

    vec3 lab = linearBT2020ToOklab(rgb);

    float L = lab.x;
    float eps = 0.00001;
    float C = max(eps, sqrt(lab.y * lab.y + lab.z * lab.z));
    float a_ = lab.y / C;
    float b_ = lab.z / C;

    float L0 = 0.5;

    float t = find_gamut_intersection(a_, b_, L, C, L0);
    float L_clipped = L0 * (1.0 - t) + t * L;
    float C_clipped = t * C;

    return oklabToLinearDisplayP3(vec3(L_clipped, C_clipped * a_, C_clipped * b_));
}

vec3 gamut_clip_project_to_L_cusp(vec3 rgb) {
    if (rgb.r < 1.0 && rgb.g < 1.0 && rgb.b < 1.0 && rgb.r > 0.0 && rgb.g > 0.0 && rgb.b > 0.0)
    return rgb;

    vec3 lab = linearBT2020ToOklab(rgb);

    float L = lab.x;
    float eps = 0.00001;
    float C = max(eps, sqrt(lab.y * lab.y + lab.z * lab.z));
    float a_ = lab.y / C;
    float b_ = lab.z / C;

    // The cusp is computed here and in find_gamut_intersection, an optimized solution would only compute it once.
    vec2 cusp = find_cusp(a_, b_);
    float cuspL = cusp.x;
    float cuspC = cusp.y;

    float L0 = cuspL;

    float t = find_gamut_intersection(a_, b_, L, C, L0);

    float L_clipped = L0 * (1.0 - t) + t * L;
    float C_clipped = t * C;

    return oklabToLinearDisplayP3(vec3(L_clipped, C_clipped * a_, C_clipped * b_));
}

vec3 gamut_clip_adaptive_L0_0_5(vec3 rgb, float alpha) {
    if (rgb.r < 1.0 && rgb.g < 1.0 && rgb.b < 1.0 && rgb.r > 0.0 && rgb.g > 0.0 && rgb.b > 0.0)
    return rgb;

    vec3 lab = linearBT2020ToOklab(rgb);

    float L = lab.x;
    float eps = 0.00001;
    float C = max(eps, sqrt(lab.y * lab.y + lab.z * lab.z));
    float a_ = lab.y / C;
    float b_ = lab.z / C;

    float Ld = L - 0.5;
    float e1 = 0.5 + abs(Ld) + alpha * C;
    float L0 = 0.5 * (1.0 + sign(Ld) * (e1 - sqrt(e1 * e1 - 2.0 * abs(Ld))));

    float t = find_gamut_intersection(a_, b_, L, C, L0);
    float L_clipped = L0 * (1.0 - t) + t * L;
    float C_clipped = t * C;

    return oklabToLinearDisplayP3(vec3(L_clipped, C_clipped * a_, C_clipped * b_));
}

vec3 gamut_clip_adaptive_L0_L_cusp(vec3 rgb, float alpha) {
    if (rgb.r < 1.0 && rgb.g < 1.0 && rgb.b < 1.0 && rgb.r > 0.0 && rgb.g > 0.0 && rgb.b > 0.0)
    return rgb;

    vec3 lab = linearBT2020ToOklab(rgb);

    float L = lab.x;
    float eps = 0.00001;
    float C = max(eps, sqrt(lab.y * lab.y + lab.z * lab.z));
    float a_ = lab.y / C;
    float b_ = lab.z / C;

    // The cusp is computed here and in find_gamut_intersection, an optimized solution would only compute it once.
    vec2 cusp = find_cusp(a_, b_);
    float cuspL = cusp.x;
    float cuspC = cusp.y;

    float Ld = L - cuspL;
    float k = 2.0 * (Ld > 0.0 ? 1.0 - cuspL : cuspL);

    float e1 = 0.5 * k + abs(Ld) + alpha * C / k;
    float L0 = cuspL + 0.5 * (sign(Ld) * (e1 - sqrt(e1 * e1 - 2.0 * k * abs(Ld))));

    float t = find_gamut_intersection(a_, b_, L, C, L0);
    float L_clipped = L0 * (1.0 - t) + t * L;
    float C_clipped = t * C;

    return oklabToLinearDisplayP3(vec3(L_clipped, C_clipped * a_, C_clipped * b_));
}

mat3 uInputTransformMatrix = M_BT2020_TO_XYZ;
mat3 uOutputTransformMatrix = M_XYZ_TO_BT2020;

//
// Tone mapping
//
highp vec3 InputTransform(const highp vec3 color) {
    return uInputTransformMatrix * color;
}

highp vec3 OutputTransform(const highp vec3 color) {
    return uOutputTransformMatrix * color;
}

const float displayMaxLuminance = 800.0;
const float maxMasteringLuminance = 1000.0;
const float maxContentLuminance = 1000.0;

// Perceptual tone mapping curve (EETF) specified in ITU-R Report BT.2390.
// This is the recommended curve to use for typical HDR-mastered content.
highp vec3 ToneMap(highp vec3 color, float maxInLumi, float maxOutLumi) {
//    float maxMasteringLumi = maxMasteringLuminance;
//    float maxContentLumi = maxContentLuminance;
//    float maxInLumi = min(maxMasteringLumi, maxContentLumi);
//    float maxOutLumi = displayMaxLuminance;

//    highp float nits = color.r * 0.2627 + color.g * 0.6780 + color.b * 0.0593;
    float nits = color.y;

    // clamp to max input luminance
    nits = clamp(nits, 0.0, maxInLumi);

    // scale [0.0, maxInLumi] to [0.0, maxOutLumi]
    if (maxInLumi <= maxOutLumi) {
        return color * (maxOutLumi / maxInLumi);
    } else {
        // three control points
        const float x0 = 10.0;
        const float y0 = 17.0;
        float x1 = maxOutLumi * 0.75;
        float y1 = x1;
        float x2 = x1 + (maxInLumi - x1) / 2.0;
        float y2 = y1 + (maxOutLumi - y1) * 0.75;

        // horizontal distances between the last three control points
        float h12 = x2 - x1;
        float h23 = maxInLumi - x2;
        // tangents at the last three control points
        float m1 = (y2 - y1) / h12;
        float m3 = (maxOutLumi - y2) / h23;
        float m2 = (m1 + m3) / 2.0;

        if (nits < x0) {
            // scale [0.0, x0] to [0.0, y0] linearly
            float slope = y0 / x0;
            return color * slope;
        } else if (nits < x1) {
            // scale [x0, x1] to [y0, y1] linearly
            float slope = (y1 - y0) / (x1 - x0);
            nits = y0 + (nits - x0) * slope;
        } else if (nits < x2) {
            // scale [x1, x2] to [y1, y2] using Hermite interp
            float t = (nits - x1) / h12;
            nits = (y1 * (1.0 + 2.0 * t) + h12 * m1 * t) * (1.0 - t) * (1.0 - t) +
            (y2 * (3.0 - 2.0 * t) + h12 * m2 * (t - 1.0)) * t * t;
        } else {
            // scale [x2, maxInLumi] to [y2, maxOutLumi] using Hermite interp
            float t = (nits - x2) / h23;
            nits = (y2 * (1.0 + 2.0 * t) + h23 * m2 * t) * (1.0 - t) * (1.0 - t) +
            (maxOutLumi * (3.0 - 2.0 * t) + h23 * m3 * (t - 1.0)) * t * t;
        }
    }

    // color.y is greater than x0 and is thus non-zero
    return color * (nits / color.y);
}

highp vec3 ScaleLuminance(highp vec3 color) {
    // The formula is:
    // alpha * pow(Y, lambada - 1.0) * color + beta;
    // where alpha is 1000.0, lambada is 1.2, beta is 0.0.
    return color * 1000.0 * pow(color.y, 0.2);
}

highp vec3 NormalizeLuminance(highp vec3 color) {
    return color / displayMaxLuminance;
}

highp vec3 OOTF(const highp vec3 color) {
    return NormalizeLuminance(ToneMap(ScaleLuminance(color), 1000.0, 800.0));
}

float apply_bt2390(float x, const float maxLum) {
    float ks = 1.5 * maxLum - 0.5;
    float tb = (x - ks) / (1.0 - ks);
    float tb2 = tb * tb;
    float tb3 = tb2 * tb;
    float pb = (2.0 * tb3 - 3.0 * tb2 + 1.0) * ks
    + (tb3 - 2.0 * tb2 + tb) * (1.0 - ks)
    + (-2.0 * tb3 + 3.0 * tb2) * maxLum;
    //x = mix(pb, x, lessThan(x, ks));
    x = (x < ks) ? x : pb;
    return x;
}

//
// BT.2446 A
//
highp vec3 BT2446_A(highp vec3 rgb, int hdrMode, float Lhdr, float Lsdr) {
    bool PQMode = hdrMode == 1;
    if (!PQMode) {
        float lambada = 1.2 + 0.42 * log10(Lhdr / 1000.0);
        rgb = HLG_OOTF_3f(EOTF(rgb), lambada);
    } else {
        vec3 xyz = rgb;
//        xyz.x = xyz.x * 1000.0; // dst peak
//        xyz.y = xyz.y * 1000.0;
//        xyz.z = xyz.z * 1000.0;
//
//        float sig_peak_pq = Y_2_ST2084(1500.0);
//        float scale = 1.0 / sig_peak_pq;
//
//        xyz.x = Y_2_ST2084(xyz.x) * scale;
//        xyz.y = Y_2_ST2084(xyz.y) * scale;
//        xyz.z = Y_2_ST2084(xyz.z) * scale;
//
//        float maxLum = Y_2_ST2084(1000.0) * scale;
//        xyz.x = apply_bt2390(xyz.x, maxLum) * sig_peak_pq;
//        xyz.y = apply_bt2390(xyz.y, maxLum) * sig_peak_pq;
//        xyz.z = apply_bt2390(xyz.z, maxLum) * sig_peak_pq;

        xyz.x = ST2084_2_Y(xyz.x) * pq_C / Lhdr;
        xyz.y = ST2084_2_Y(xyz.y) * pq_C / Lhdr;
        xyz.z = ST2084_2_Y(xyz.z) * pq_C / Lhdr;

        xyz.x = xyz.x / pq_C;
        xyz.y = xyz.y / pq_C;
        xyz.z = xyz.z / pq_C;

        rgb = xyz;
        Lhdr = 1000.0;
    }

    float CoeffAdj=1.0;
    float alpha=1.0/2.404;
    float pHDR=1.0+32.0*pow(Lhdr/10000.0,alpha);
    float ilogpHDR=1.0/log(pHDR);
    float pSDR=1.0+32.0*pow(Lsdr/10000.0,alpha);
    float ipSDR=1.0/(pSDR-1.0);
    float Coeffb=1.0/1.8814;
    float Coeffr=1.0/1.4746;

    highp float R = min(pow(CoeffAdj*rgb.r,alpha),1.0);
    highp float G = min(pow(CoeffAdj*rgb.g,alpha),1.0);
    highp float B = min(pow(CoeffAdj*rgb.b,alpha),1.0);
    highp float Y = R * 0.2627 + G * 0.6780 + B * 0.0593;

    highp float y=ilogpHDR*log(1.0+(pHDR-1.0)*Y);
    if (y<=0.7399) y*=1.0770;
    else
    {
        if (y<0.9909) y=(2.7811-1.1510*y)*y-0.6302;
        else y=0.5000*y+0.5000;
    }
    y=ipSDR*(pow(pSDR,y)-1.0);

    highp float f=(Y==0.0)?0.0:y/(1.1*Y);

    highp float Cb=f*Coeffb*(B-Y);
    highp float Cr=f*Coeffr*(R-Y);
    highp float Yt=y-max(0.1*Cr, 0.0);

    const highp float a = 0.2627;
    const highp float b = 0.6780;
    const highp float c = 0.0593;  // a + b + c = 1
    const highp float d = 1.8814;  // 2 * (a + b)
    const highp float e = 1.4747;  // 2 * (1 - a)

    rgb.r = Yt + e * Cr;
    rgb.g = Yt - (a * e / b) * Cr - (c * d / b) * Cb;
    rgb.b = Yt + d * Cb;

    vec3 XYZ = M_BT2020_TO_XYZ * rgb;
    rgb = M_XYZ_TO_DISPLAY_P3 * XYZ;

    rgb.r = pow(rgb.r, 2.404);
    rgb.g = pow(rgb.g, 2.404);
    rgb.b = pow(rgb.b, 2.404);

    return clamp(rgb, 0.0, 1.0);
}

//
// IRU BT.2446 C
//

highp vec3 BT2446_C(highp vec3 rgb, int hdrMode, float Lhdr, float Lsdr) {
    bool PQMode = hdrMode == 1;

    const bool ChromaC = false;

    float hdrScale = (PQMode) ? pq_C : Lhdr;

    float pct_ref = (PQMode) ? 0.58 : 0.75;
    float pct_ip = 0.80;
    float pct_wp = 0.96;

    float pct_sdr_skin = 0.70;
    float pct_hdr_skin = (PQMode) ? 0.44 : 0.50;

    //  alpha=0.15, sigma=0.5

    // parameter of chrosstalk matrix used to reduce the chroma
    // before applying the tonecurve.
    float alpha = 0.15;

    // parameter of chroma reduction.
    float sigma = 0.5;

    // ===================================================
    // 1. Conversion to linear display light signals
    // ===================================================

    vec3 RGB = rgb;

    float lambada;
    if (!PQMode) {
        lambada = 1.2 + 0.42 * log10(Lhdr / 1000.0);
        RGB = HLG_OOTF_3f(EOTF(RGB), lambada);
    } else {
        lambada = 1.0;
        RGB = ST2084_2_Y_f3(RGB) / pq_C;
    }

    // ===================================================
    // 2. Crosstalk matrix
    // ===================================================

    RGB = crosstalk(RGB, alpha);

    // ===================================================
    // 3. Conversion to XYZ
    // ===================================================

    vec3 XYZ = InputTransform(RGB);

    // ===================================================
    // 4 Calculate k1 ... k4 parameters
    // ===================================================

    // D65 white point
    const float Wx = 0.3127;
    const float Wy = 0.3290;

    float k1;
    if (PQMode) {
        k1 = Lsdr * pow(pct_sdr_skin, 2.404) / (Lhdr * ST2084_2_Y(pct_hdr_skin) / pq_C);
    } else {
        k1 = Lsdr * pow(pct_sdr_skin, 2.404) / (Lhdr * HLG_OOTF(HLG_inv_OETF(pct_hdr_skin), lambada));
    }

    float Ysdr_ip, Ysdr_wp;
    Ysdr_ip = Lsdr * pow(pct_ip, 2.404);
    Ysdr_wp = Lsdr * pow(pct_wp, 2.404);

    float Yhdr_ip, Yhdr_ref;
    if (PQMode) {
        Yhdr_ref = Lhdr * ST2084_2_Y(pct_ref) / pq_C;
    } else {
        Yhdr_ref = Lhdr * HLG_OOTF(HLG_inv_OETF(pct_ref), lambada);
    }

    Yhdr_ip = Ysdr_ip / k1;

    // How to calculate k3?
    float k3 = 0.7409841362;
    float k2 = k1 * (1.0 - k3) * Yhdr_ip;
    float k4 = k1 * Yhdr_ip - k2 * log(1.0 - k3);

    const float Xmin = 0.0000000000;
    const float Xmax = 3127.0 / 3290.0; // 0.9504284859;

    const float Ymin = 0.0000000000;
    const float Ymax = 1.0000000000;

    const float Zmin = 0.0000000000;
    const float Zmax = 3583.0 / 3290.0; // 1.0889004469;

    const float CoeffX = Xmax - Xmin;
    const float CoeffY = Ymax - Ymin;
    const float CoeffZ = Zmax - Zmin;

    const float iCoeffX = 1.0 / CoeffX;
    const float iCoeffY = 1.0 / CoeffY;
    const float iCoeffZ = 1.0 / CoeffZ;

    // Xn: 193.0815944167
    float Xn = Yhdr_ref * Wx / Wy;
    // Yn: 15.1479794833
    float Yn = Yhdr_ref;
    // Zn: 221.2124632930
    float Zn = Yhdr_ref * (1.0 - Wx - Wy) / Wy;

//    float x = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
//    float y = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);
//    float z = 1.0 - x -y;
//    float Y = XYZ.y * hdrScale;

    float x = XYZ.x;// * hdrScale;
    float y = XYZ.y;// * hdrScale;
    float z = XYZ.z;// * hdrScale;

    // ===================================================
    // 4. Conversion of XYZ to CIE L*a*b*
    // ===================================================

    float L, a, b, fx, fy, fz, t;

    float d1 = pow(6.0 / 29.0, 3.0);
    float d2 = pow(29.0 / 6.0, 2.0) / 3.0;
    float d4 = 4.0 / 29.0;
    float pcoeff = 1.0 / 3.0;

    t = x * hdrScale / Xn;
    fx = (t > d1) ? pow(t, pcoeff) : d2 * t + d4;

    t = y * hdrScale / Yn;
    fy = (t > d1) ? pow(t, pcoeff) : d2 * t + d4;

    t = z *hdrScale / Zn;
    fz = (t > d1) ? pow(t, pcoeff) : d2 * t + d4;

    L = 116.0 * fy - 16.0;
    a = 500.0 * (fx - fy);
    b = 200.0 * (fy - fz);

    // ===================================================
    // 5. Chroma correction above HDR reference White
    // ===================================================

    float Cab, hab, fcor;

    float Lref = 116.0 - 16.0;
    float Lmax = 116.0 * pow(Lhdr / Yn, pcoeff) - 16.0; // Lhdr / Yn => 1000.0+ / 203?
    float iLmr = 1.0 / (Lmax - Lref);

    Cab = sqrt(a * a + b * b);
    hab = atan2(b, a);

    fcor = (L > Lref) ? (1.0 - sigma * (L - Lref) * iLmr) : 1.0;
    // must not be negative values
    fcor = max(fcor, 0.0);
    Cab *= fcor;
    a = Cab * cos(hab);
    b = Cab * sin(hab);

    // ===================================================
    // 6. Conversion of CIE L*a*b* to XYZ
    // ===================================================

    float sXn = Xn / hdrScale;
    float sYn = Yn / hdrScale;
    float sZn = Zn / hdrScale;

    float i116 = 1.0 / 116.0;
    float i500 = 1.0 / 500.0;
    float i200 = 1.0 / 200.0;

    float delta = 6.0 / 29.0;
    float d3 = 3.0 * pow(6.0 / 29.0, 2.0);

    fy = (L + 16.0) * i116;
    fx = fy + a * i500;
    fz = fy - b * i200;

    x = (fx > delta) ? sXn * pow(fx, 3.0) : (fx - d4) * d3 * sXn;
    //x = clamp((x - Xmin) * iCoeffX, 0.0, 1.0);

    y = (fy > delta) ? sYn * pow(fy, 3.0) : (fy - d4) * d3 * sYn;
    //y = clamp((y - Ymin) * iCoeffY, 0.0, 1.0);

    z = (fz > delta) ? sZn * pow(fz, 3.0) : (fz - d4) * d3 * sZn;
   // z = clamp((z - Zmin) * iCoeffZ, 0.0, 1.0);

    // ===================================================
    // 7. Tone mapping
    // ===================================================

    float yhdr = y * hdrScale;
    float ysdr;
    if (yhdr < Yhdr_ip /* 69 */) {
        ysdr = k1 * yhdr / Lsdr;
    } else {
        ysdr = (k4 + k2 * log(yhdr / Yhdr_ip - k3)) / Lsdr;
    }

    float iy = (y == 0.0) ? 0.0 : ysdr / y;

    y = ysdr;
    x *= iy;
    z *= iy;


//    // ===================================================
//    // 7. Tone mapping
//    // ===================================================
//
//    float k1 = 0.83802;
//    float k2 = 15.09968;
//    float k3 = 0.74204;
//    float k4 = 78.99439;
//    float Yhdr_ip = 58.5 / k1;
//    float yhdr = Y;
//    float ysdr;
//    if (yhdr < Yhdr_ip /* 69 */) {
//        ysdr = k1 * yhdr;
//    } else {
//        ysdr = k2 * log(yhdr / Yhdr_ip - k3) + k4;
//    }
//
//    float iy = (y == 0.0) ? 0.0 : ysdr / y;
//
//    y = ysdr;
//    x *= iy;
//    z *= iy;

    XYZ = clamp(vec3(x, y, z) , 0.0, 1.0);

    // ===================================================
    // 8. Color gamut mapping
    // ===================================================

    RGB = M_XYZ_TO_DISPLAY_P3 * XYZ;
//    RGB = inverse_crosstalk(RGB, alpha);


//     RGB = gamut_clip_preserve_lightness(RGB);
//     RGB = gamut_clip_project_to_0_5(RGB);
//     RGB = gamut_clip_project_to_L_cusp(RGB);
//    RGB = gamut_clip_adaptive_L0_0_5(RGB, 0.05);
//     RGB = gamut_clip_adaptive_L0_L_cusp(RGB, 0.05);
//    RGB.r = pow(RGB.r, 1.0/2.404);
//    RGB.g = pow(RGB.g, 1.0/2.404);
//    RGB.b = pow(RGB.b, 1.0/2.404);
//    RGB = OETF_sRGB(RGB);
    return clamp(RGB, 0.0, 1.0);
}

highp vec3 BT2446_C_Simple(highp vec3 rgb, int hdrMode, float Lhdr, float Lsdr) {
    bool PQMode = hdrMode == 1;

    float alpha = 0.0;
    float hdrScale = (PQMode) ? pq_C : Lhdr;

    vec3 RGB = crosstalk(rgb, alpha);

    float lambada;
    if (!PQMode) {
        lambada = 1.2 + 0.42 * log10(Lhdr / 1000.0);
        RGB = HLG_OOTF_3f(EOTF(RGB), lambada);
    } else {
        lambada = 1.0;
        RGB = ST2084_2_Y_f3(RGB) / pq_C;
    }

    vec3 XYZ = InputTransform(RGB);

    float x = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
    float y = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);
    float z = 1.0 - x -y;
    float Y = XYZ.y * hdrScale;

    float k1 = 0.83802;
    float k2 = 15.09968;
    float k3 = 0.74204;
    float k4 = 78.99439;
    float Yhdr_ip = 58.5 / k1;
    float yhdr = Y;
    float ysdr;
    if (yhdr < Yhdr_ip /* 69 */) {
        ysdr = k1 * yhdr;
    } else {
        ysdr = k2 * log(yhdr / Yhdr_ip - k3) + k4;
    }

    float iy = (y == 0.0) ? 0.0 : ysdr / y;

    y = ysdr;
    x *= iy;
    z *= iy;

    XYZ = clamp(vec3(x, y, z) / Lsdr, 0.0, 1.0);


//    RGB = M_XYZ_TO_BT2020 * XYZ;
//    RGB = inverse_crosstalk(RGB, alpha);
//    RGB = gamut_clip_preserve_lightness(RGB);
//    RGB = gamut_clip_project_to_0_5(RGB);
//    RGB = gamut_clip_project_to_L_cusp(RGB);
//    RGB = gamut_clip_adaptive_L0_0_5(RGB, 0.05);
//    RGB = gamut_clip_adaptive_L0_L_cusp(RGB, 0.06);
//    RGB = gamut_clip_adaptive_L0_0_5(RGB, 0.5);
    RGB = M_XYZ_TO_DISPLAY_P3 * XYZ;

    //    RGB.r = pow(RGB.r, 1.0/2.404);
    //    RGB.g = pow(RGB.g, 1.0/2.404);
    //    RGB.b = pow(RGB.b, 1.0/2.404);
    //    RGB = OETF_sRGB(RGB);
    return clamp(RGB, 0.0, 1.0);
}


// HLG
//const int uHdrMode = 2;
//float uHdrMaxMasteringLuminance = 1000.0;
//float uSdrWhiteTargetingLuminance = 100.0;

// PQ
const int uHdrMode = 1;
float uHdrMaxMasteringLuminance = 10000.0;
float uSdrWhiteTargetingLuminance = 200.0;

void main() {
    gl_FragColor = texture2D(sTexture, vTextureCoord);

//    if (0.0 <= gl_FragColor.r && gl_FragColor.r <= 1.0) {
//        gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);
//    } else {
//        gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
//    }

    if (uHdrMode > 0) {
        // method 1
        // gl_FragColor.rgb = OETF(OutputTransform(InputTransform(EOTF(gl_FragColor.rgb))));

        // method 2
        // gl_FragColor.rgb = OETF_sRGB(OutputTransform(OOTF(InputTransform(EOTF(gl_FragColor.rgb)))));

        // BT2446_C
//         gl_FragColor.rgb = BT2446_C(gl_FragColor.rgb, uHdrMode, uHdrMaxMasteringLuminance, uSdrWhiteTargetingLuminance);
        gl_FragColor.rgb = BT2446_C_Simple(gl_FragColor.rgb, uHdrMode, uHdrMaxMasteringLuminance, uSdrWhiteTargetingLuminance);

        // BT2446_A
//         gl_FragColor.rgb = BT2446_A(gl_FragColor.rgb, uHdrMode, uHdrMaxMasteringLuminance, uSdrWhiteTargetingLuminance);
        gl_FragColor.a = 1.0;
    }
}