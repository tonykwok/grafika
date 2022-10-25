package com.android.grafika.gles;

import android.support.annotation.IntDef;
import android.support.annotation.NonNull;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class ColorSpace {
    @IntDef({TEX_GL_COLORSPACE_SRGB, TEX_GL_COLORSPACE_DISPLAY_P3,
            TEX_GL_COLORSPACE_BT2020_HLG, TEX_GL_COLORSPACE_BT2020_PQ})
    @Retention(RetentionPolicy.SOURCE)
    public @interface TexColorSpace {}
    public static final int TEX_GL_COLORSPACE_SRGB = 0;
    public static final int TEX_GL_COLORSPACE_DISPLAY_P3 = 1;
    public static final int TEX_GL_COLORSPACE_BT2020_HLG = 2;
    public static final int TEX_GL_COLORSPACE_BT2020_PQ = 3;

    @IntDef({EGL_GL_COLORSPACE_LINEAR_KHR, EGL_GL_COLORSPACE_SRGB_KHR,
            EGL_GL_COLORSPACE_DISPLAY_P3_EXT, EGL_GL_COLORSPACE_DISPLAY_P3_LINEAR_EXT,
            EGL_GL_COLORSPACE_DISPLAY_P3_PASSTHROUGH_EXT, EGL_GL_COLORSPACE_BT2020_PQ_EXT,
            EGL_GL_COLORSPACE_BT2020_LINEAR_EXT})
    @Retention(RetentionPolicy.SOURCE)
    public @interface EglColorSpace {}
    public static final int EGL_GL_COLORSPACE_LINEAR_KHR = 0x308A;
    // https://www.khronos.org/registry/EGL/extensions/KHR/EGL_KHR_gl_colorspace.txt
    public static final int EGL_GL_COLORSPACE_SRGB_KHR = 0x3089;
    // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_display_p3.txt
    public static final int EGL_GL_COLORSPACE_DISPLAY_P3_EXT = 0x3363;
    public static final int EGL_GL_COLORSPACE_DISPLAY_P3_LINEAR_EXT = 0x3362;
    // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_display_p3_passthrough.txt
    public static final int EGL_GL_COLORSPACE_DISPLAY_P3_PASSTHROUGH_EXT = 0x3490;
    // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_bt2020_linear.txt
    public static final int EGL_GL_COLORSPACE_BT2020_PQ_EXT = 0x3340;
    public static final int EGL_GL_COLORSPACE_BT2020_LINEAR_EXT = 0x333F;

    public static Set<String> getEglColorSpaceExtensions(@EglColorSpace int eglColorSpace) {
        final Set<String> requiredExtensions = new HashSet<>();
        requiredExtensions.add("EGL_KHR_gl_colorspace");
        switch (eglColorSpace) {
            // The default value of EGL_GL_COLORSPACE_KHR is EGL_GL_COLORSPACE_LINEAR_KHR.
            case EGL_GL_COLORSPACE_LINEAR_KHR:
                break;
            // https://www.khronos.org/registry/EGL/extensions/KHR/EGL_KHR_gl_colorspace.txt
            case EGL_GL_COLORSPACE_SRGB_KHR:
                break;
            // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_display_p3.txt
            case EGL_GL_COLORSPACE_DISPLAY_P3_EXT:
                requiredExtensions.add("EGL_EXT_gl_colorspace_display_p3");
                break;
            case EGL_GL_COLORSPACE_DISPLAY_P3_LINEAR_EXT:
                requiredExtensions.add("EGL_EXT_gl_colorspace_display_p3_linear");
                break;
            // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_display_p3_passthrough.txt
            case EGL_GL_COLORSPACE_DISPLAY_P3_PASSTHROUGH_EXT:
                requiredExtensions.add("EGL_EXT_gl_colorspace_display_p3_passthrough");
                break;
            // https://www.khronos.org/registry/EGL/extensions/EXT/EGL_EXT_gl_colorspace_bt2020_linear.txt
            case EGL_GL_COLORSPACE_BT2020_PQ_EXT:
                requiredExtensions.add("EGL_EXT_gl_colorspace_bt2020_pq");
                break;
            case EGL_GL_COLORSPACE_BT2020_LINEAR_EXT:
                requiredExtensions.add("EGL_EXT_gl_colorspace_bt2020_linear");
                break;
            default:
                throw new IllegalArgumentException("Unknown EglColorSpace: " + eglColorSpace);
        }
        return Collections.unmodifiableSet(requiredExtensions);
    }

    public final @TexColorSpace int texColorSpace;
    public final @EglColorSpace int eglColorSpace;

    public static final ColorSpace sDefaultColorSpace = new ColorSpace(
            TEX_GL_COLORSPACE_SRGB,
            EGL_GL_COLORSPACE_LINEAR_KHR
    );

    public ColorSpace(@TexColorSpace int texColorSpace, @EglColorSpace int eglColorSpace) {
        this.texColorSpace = texColorSpace;
        this.eglColorSpace = eglColorSpace;
    }

    @Override
    public @NonNull String toString() {
        return "ColorSpace(tex: " + texColorSpace + ", egl: " + eglColorSpace + ")";
    }
}
