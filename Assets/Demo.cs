using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Demo : MonoBehaviour {
	Texture2D img;

	MSAFluidSolver fluidSolver;
	MSAFluidDrawer fluidDrawer;

	Texture2D texture;
	Texture2D brushTexture;
	Vector2 lastHitPos;

	 float                   colorMult = 2f;
	 float                   velocityMult = 20;
	 int                     fluidCellsX = 55;
	 bool                    resizeFluid;
	 bool                    drawParticles = true;
	bool                    drawFluid = true;

	Color currColor = Color.white;
	float hue;
	public ParticleSystem ps;
	public Camera renderCamera;

	float fpsTime;
	int fpsFrame;
	string lastFps;

	void Start () {

		fluidSolver = new MSAFluidSolver();
		fluidDrawer = new MSAFluidDrawer();

		colorMult = 2f;
		fluidCellsX = 80;

		fluidSolver.setup(fluidCellsX,fluidCellsX);
		fluidSolver.enableRGB(true).setFadeSpeed(0.003f).setDeltaT(0.5f).setVisc(0.0001f);
		fluidDrawer.setup(fluidSolver);

		texture = new Texture2D(fluidCellsX,fluidCellsX);
		GetComponent<Renderer>().material.mainTexture = texture;
		
		currColor = ColorFromHSV(Rand.randFloat() * 360f,0.75f,0.75f);
		ps.startColor = currColor;
	}
	
	void mouseMoved() {   
		RaycastHit hit;
		Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
		
		if(Physics.Raycast(ray, out hit))
		{
			Vector2 hitPos = hit.textureCoord;
			Vector2 mouseVel = hitPos - lastHitPos;
			
			if(hit.collider == GetComponent<Collider>())
			{
				addToFluid(hit.point , hitPos, mouseVel, true,true);
				lastHitPos = hitPos;
			}
		}
	}
    
	void mouseDragged() {   
    	RaycastHit hit;
        Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
        
        if(Physics.Raycast(ray, out hit))
		{
			Vector2 hitPos = hit.textureCoord;
			Vector2 mouseVel = hitPos - lastHitPos;

			if(hit.collider == GetComponent<Collider>())
			{
				addToFluid(hit.point , hitPos, mouseVel, false,true);
				lastHitPos = hitPos;
			}
        } 
   }
    

	void Update () {
		if(resizeFluid) 	{
			fluidSolver.setSize(fluidCellsX, fluidCellsX);
			fluidDrawer.setup(fluidSolver);

			Destroy(texture);
			texture = new Texture2D(fluidCellsX,fluidCellsX);
			GetComponent<Renderer>().material.mainTexture = texture;

			resizeFluid = false;
		}
		
		UpdatePartical();
		fluidSolver.update();

		if(drawFluid)
		{
			fluidDrawer.draw();

			Color[] cols = texture.GetPixels();
			for( int i = 0; i < cols.Length; ++i ) {
				cols[i] = fluidDrawer._pixels[i];
			}

			texture.SetPixels(cols);
			texture.Apply();
		}	

		if(Input.GetMouseButton(0))
			mouseDragged();
		else 
			mouseMoved();

		if(Input.GetMouseButtonDown(0))
			hue = Rand.randFloat() * 360f;
		else 
		{
			hue += Time.deltaTime * 18f;
			if(hue >= 360f)
				hue = 0;		
		}
		currColor = ColorFromHSV(hue,0.8f,0.8f);
		ps.startColor = currColor;

		fpsTime += Time.deltaTime;
		if(fpsTime >= 1f)
		{
			lastFps = ((float)fpsFrame / fpsTime).ToString("##.##");
			fpsTime = 0f;
			fpsFrame = 0;
		}
		fpsFrame++;
	}

	public static Color ColorFromHSV(float hue, float saturation, float value)
	{
		int hi = Mathf.FloorToInt(hue / 60) % 6;
		double f = hue / 60 - Mathf.FloorToInt(hue / 60);
		
		value = value * 255;
		int v = (int)(value);
		int p = (int)(value * (1 - saturation));
		int q = (int)(value * (1 - f * saturation));
		int t = (int)(value * (1 - (1 - f) * saturation));
		
		if (hi == 0)
			return new Color( v / 255f, t/255f, p/255f);
		else if (hi == 1)
			return new Color(q / 255f, v / 255f, p / 255f);
		else if (hi == 2)
			return new Color( p / 255f, v / 255f, t / 255f);
		else if (hi == 3)
			return new Color(p / 255f, q / 255f, v / 255f);
		else if (hi == 4)
			return new Color(t / 255f, p / 255f, v / 255f);
		else
			return new Color(v / 255f, p / 255f, q / 255f);
	}

	void OnGUI()
	{
		GUI.backgroundColor = Color.gray;
		GUILayout.BeginVertical();
		addSlider("fluidCellsX", ref fluidCellsX, 20, 400);
		addButton("resizeFluid",ref  resizeFluid);
		addSlider("colorMult",ref  colorMult, 0f, 2.5f);
		addSlider("velocityMult",ref  velocityMult, 0, 100);
		addSlider("fs.viscocity",ref  fluidSolver.viscocity, 0.0f, 0.01f);
		addSlider("fs.colorDiffusion",ref  fluidSolver.colorDiffusion, 0.0f, 0.0003f); 
		addSlider("fs.fadeSpeed",ref  fluidSolver.fadeSpeed, 0.0f, 0.1f); 
		addSlider("fs.solverIterations",ref  fluidSolver.solverIterations, 1, 50); 
		addSlider("fs.deltaT",ref  fluidSolver.deltaT, 0.1f, 5f);
		addComboBox("fd.drawMode", ref fluidDrawer.drawMode);
		addToggle("fs.doRGB", ref fluidSolver.doRGB); 
		addToggle("fs.doVorticityConfinement",ref fluidSolver.doVorticityConfinement); 
		addToggle("drawFluid",ref drawFluid); 
		addToggle("drawParticles",ref drawParticles); 
		addToggle("fs.wrapX",ref fluidSolver.wrap_x);
		addToggle("fs.wrapY", ref fluidSolver.wrap_y);
		GUILayout.Label("FPS: " + lastFps);
		GUILayout.EndVertical();
	}

	void addSlider (string title, ref int data, int min, int max)
	{
		GUILayout.BeginHorizontal();
		GUILayout.Label(title);
		data = (int)GUILayout.HorizontalSlider(data,min,max,GUILayout.Width(100f));
		GUILayout.EndHorizontal();
		GUILayout.BeginHorizontal();
		GUILayout.FlexibleSpace();
		GUILayout.Label(" (" + data + ")");
		GUILayout.EndHorizontal();
	}

	void addSlider (string title, ref float data, float min, float max)
	{
		GUILayout.BeginHorizontal();
		GUILayout.Label(title,GUILayout.Width(100f));
		data = GUILayout.HorizontalSlider(data,min,max,GUILayout.Width(130f));

		GUILayout.Label(" (" + data + ")");
		GUILayout.EndHorizontal();
	}

	void addButton (string title, ref bool data)
	{
		GUILayout.BeginHorizontal();
		GUILayout.Label(title);
		data = GUILayout.Button(title);
		GUILayout.EndHorizontal();
	}

	void addComboBox (string title, ref MSAFluidDrawer.DrawMode drawMode)
	{
		GUILayout.BeginHorizontal();
		GUILayout.Label(title);
		string current = drawMode.ToString();
		bool c  = GUILayout.Button(drawMode.ToString());
		if(c)
		{
			int i = (int)drawMode;
			i++;
			if(i >= (int	)MSAFluidDrawer.DrawMode.kDrawCount)
				i = 0;
			drawMode = (MSAFluidDrawer.DrawMode)i;
		}
		GUILayout.EndHorizontal();
	}

	void addToggle (string title, ref bool doRGB)
	{
		GUILayout.BeginHorizontal();
		GUILayout.Label(title);
		doRGB = GUILayout.Toggle(doRGB,title);
		GUILayout.EndHorizontal();
	}

	void addToFluid(Vector3 worldPos ,Vector2 pos,Vector2 vel,bool addColor, bool addForce) {
		float speed = vel.x * vel.x  + vel.y * vel.y;    // balance the x and y components of speed with the screen aspect ratio
		if(speed == 0)
			return;
		vel *= 2f;

		int index = fluidSolver.getIndexForNormalizedPosition(pos.x, pos.y);

		if(drawFluid)
		{
			if(addColor)
			{
				Vector3 col = new Vector3(currColor.r,currColor.g,currColor.b);
				fluidSolver.addColorAtIndex(index,col * colorMult);		
			}
		}

		fluidSolver.addForceAtIndex(index,vel * velocityMult);

		if(drawParticles)
		{
			worldPos.z = 0;
			addParticles(worldPos,vel,5);
		}
	}

	void addParticles(Vector3 pos,Vector3 toDirection, int count){
		for(int i=0; i<count; i++)
		{
			pos.z = 0f;
			ps.transform.position = pos;
			ps.Emit(count);
		}
	}

	private void UpdatePartical(){

		ParticleSystem.Particle[] _particles = new ParticleSystem.Particle[ps.particleCount];
		ps.GetParticles(_particles);

		float helfsize = transform.lossyScale.x * 5f;
		float Left = -helfsize  + transform.position.x;
		float Right = helfsize + transform.position.x;
		float top = helfsize + transform.position.y;
		float bottom = -helfsize + transform.position.y;

		for (int i = 0; i < _particles.Length; i++){
			Vector3 diff = _particles[i].position ;
			float xCoord = (diff.x - Left) / (Right - Left);
			float yCoord = (diff.y - Left) / (Right - Left);

			int index = fluidSolver.getIndexForNormalizedPosition(xCoord,yCoord);

			Vector2 vel = fluidSolver.getVelocityAtIndex(index) * 10000;

			if(xCoord < 0f || xCoord > 1f)
				vel.x *= -1;
			if(yCoord < 0f || yCoord > 1f)
				vel.y *= -1;

			_particles[i].velocity = vel;
		}
		ps.SetParticles(_particles,_particles.Length);
	}
}

