varying vec2 texcoords;

uniform sampler2D texture;

void main()
{
	gl_FragColor = texture2D(texture, texcoords);
}