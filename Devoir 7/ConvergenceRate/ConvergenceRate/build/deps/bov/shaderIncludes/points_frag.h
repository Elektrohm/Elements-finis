static GLchar points_frag[]={" /*************************************************************************\n"
"  * BOV 0.1\n"
"  * A wrapper around OpenGL and GLFW (www.glfw.org) to draw simple 2D\n"
"  * graphics.\n"
"  *------------------------------------------------------------------------\n"
"  * Copyright (c) 2019-2020 C\n"
"lestin Marot <marotcelestin@gmail.com>\n"
"  *\n"
"  * This software is provided 'as-is', without any express or implied\n"
"  * warranty. In no event will the authors be held liable for any damages\n"
"  * arising from the use of this software.\n"
"  *\n"
"  * Permission is granted to anyone to use this software for any purpose,\n"
"  * including commercial applications, and to alter it and redistribute it\n"
"  * freely, subject to the following restrictions:\n"
"  *\n"
"  * 1. The origin of this software must not be misrepresented; you must not\n"
"  *    claim that you wrote the original software. If you use this software\n"
"  *    in a product, an acknowledgment in the product documentation would\n"
"  *    be appreciated but is not required.\n"
"  *\n"
"  * 2. Altered source versions must be plainly marked as such, and must not\n"
"  *    be misrepresented as being the original software.\n"
"  *\n"
"  * 3. This notice may not be removed or altered from any source\n"
"  *    distribution.\n"
"  *\n"
"  *************************************************************************/\n"
"\n"
"#version 150 core\n"
"\n"
"layout (std140) uniform objectBlock\n"
"{\n"
"	vec4 fillColor;\n"
"	vec4 outlineColor;\n"
"	vec2 localPos;\n"
"	vec2 localScale;\n"
"	float width;\n"
"	float marker;\n"
"	float outlineWidth;\n"
"	// float rotation;\n"
"	int space_type; // 0: normal sizes, 1: unzoomable, 2: unmodifable pixel size\n"
"	int colormap;\n"
"};\n"
"\n"
"in vec2 pSquare;\n"
"flat in float pixelSize;\n"
"\n"
"out vec4 outColor;\n"
"\n"
"// see https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm\n"
"\n"
"float ndot(in vec2 a, in vec2 b ) { return a.x * b.x - a.y * b.y; }\n"
"float sdRhombus( in vec2 p, in vec2 b )\n"
"{\n"
"	vec2 q = abs(p);\n"
"	float h = clamp((-2.0f * ndot(q, b) + ndot(b, b)) / dot(b, b),-1.0f, 1.0f);\n"
"	float d = length( q - 0.5f * b * vec2(1.0f - h, 1.0f + h) );\n"
"	return d * sign( q.x * b.y + q.y * b.x - b.x * b.y );\n"
"}\n"
"\n"
"// p is the coordinate, s is the size\n"
"float sdCross(in vec2 p, in vec2 s)\n"
"{\n"
"	p = abs(p); p = (p.y>p.x) \? p.yx : p.xy;\n"
"	vec2  q = p - s;\n"
"	float k = max(q.y, q.x);\n"
"	vec2  w = (k>0.0f) \? q : vec2(s.y - p.x, -k);\n"
"	return sign(k) * length(max(w, 0.0f));\n"
"}\n"
"\n"
"float sdBox(in vec2 p, in vec2 b)\n"
"{\n"
"	vec2 d = abs(p) - b;\n"
"	return length(max(d, vec2(0.0f))) + min(max(d.x, d.y), 0.0f);\n"
"}\n"
"\n"
"// signed distance to a n-star polygon with external angle en\n"
"float sdStar(in vec2 p, in float r, in float n, in float m) // m=[2,n]\n"
"{\n"
"	// these 4 lines can be precomputed for a given shape\n"
"	float an = 3.141593f / n;\n"
"	float en = 3.141593f / m;\n"
"	vec2  acs = vec2(cos(an), sin(an));\n"
"	vec2  ecs = vec2(cos(en), sin(en)); // ecs=vec2(0,1) and simplify, for regular polygon,\n"
"\n"
"	// reduce to first sector\n"
"	float bn = mod(atan(p.x, p.y), 2.0f * an) - an;\n"
"	p = length(p) * vec2(cos(bn), abs(sin(bn)));\n"
"\n"
"	// line sdf\n"
"	p -= r * acs;\n"
"	p += ecs * clamp( -dot(p, ecs), 0.0f, r * acs.y / ecs.y);\n"
"	return length(p) * sign(p.x);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"	float mwrap = mod(marker, 25.0f);      // marker%25\n"
"	float mfrac = fract(mwrap);\n"
"	float fw = mfrac*width;\n"
"\n"
"	float sdf;\n"
"	if(mwrap<1.0f)\n"
"		sdf = length(pSquare) - width;\n"
"	else if(mwrap<2.0f)\n"
"		sdf = sdBox(pSquare, vec2(fw, width));\n"
"	else if(mwrap<3.0f)\n"
"		sdf = sdBox(pSquare, vec2(width, fw));\n"
"	else if(mwrap<4.0f)\n"
"		sdf = sdBox(pSquare, vec2(fw)) + fw - width;\n"
"	else if(mwrap<5.0f)\n"
"		sdf = sdRhombus(pSquare, vec2(fw, width));\n"
"	else if(mwrap<6.0f)\n"
"		sdf = sdRhombus(pSquare, vec2(width, fw));\n"
"	else if(mwrap<7.0f)\n"
"		sdf = sdCross(pSquare, vec2(fw, 0.25f * fw)) + fw - width;\n"
"	else {\n"
"		float ad = 2.0f;   // angle divisor, between 2 and n\n"
"		float nbranch = floor(mwrap) - 4.0;\n"
"		if(mwrap<13.0f){\n"
"			ad += mfrac * mfrac * (nbranch - 2.0f);\n"
"			fw = width;\n"
"		}\n"
"		else if(mwrap<19.0f){\n"
"			nbranch -= 6.0f;\n"
"		}\n"
"		else {\n"
"			nbranch -= 12.0f;\n"
"			ad = nbranch;\n"
"		}\n"
"		sdf = sdStar(pSquare, fw, nbranch, ad) + fw - width;\n"
"	}\n"
"\n"
"	vec2 alpha = smoothstep(-pixelSize, pixelSize, -sdf - vec2(0.0f, outlineWidth));\n"
"	outColor = mix(outlineColor, fillColor, alpha.y);\n"
"	outColor.a *= alpha.x;\n"
"}\n"
};
