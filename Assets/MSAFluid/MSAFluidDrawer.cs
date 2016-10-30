using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using Vec2f = UnityEngine.Vector2;
using Vec3f = UnityEngine.Vector3;

public class MSAFluidDrawer : MSAFluidSolverUtililty
{
	public enum DrawMode{
		kDrawColor,
		kDrawMotion,
		kDrawSpeed,
		kDrawVectors,
		kDrawCount
	}
		
	static List<string> drawModeTitles;
	public List<string> getDrawModeTitles()
	{

		if(drawModeTitles.Count == 0) {
			drawModeTitles.Add("kDrawColor");
			drawModeTitles.Add("kDrawMotion");
			drawModeTitles.Add("kDrawSpeed");
			drawModeTitles.Add("kDrawVectors");
		}
		return drawModeTitles;
	}

	public bool	enabled;
	public bool	doInvert;
	public float	brightness;
	public float	velDrawThreshold;
	public float	velDrawMult;
	public int		vectorSkipCount;
			
	public DrawMode		drawMode;
			
	public Color[]		_pixels;						// pixels array to be drawn
			
	bool				_alphaEnabled;
			
	MSAFluidSolver              _fluidSolver;
	bool				_didICreateTheFluid;    // TODO: replace with shared pointer
			
	//--------------------------------------------------------------
	public MSAFluidDrawer() 
	{
		_pixels				= null;
		_fluidSolver		= null;
		_didICreateTheFluid	= false;
		
		enabled				= true;
		brightness			= 1;
		doInvert			= false;
		velDrawMult				= 1;
		vectorSkipCount		= 0;
		
		enableAlpha(false);
		
		setDrawMode(DrawMode.kDrawColor);
	}
	
	//--------------------------------------------------------------
	public MSAFluidSolver setup(int NX, int NY) 
	{
		deleteFluidSolver();
		_fluidSolver = new MSAFluidSolver();
		_fluidSolver.setup(NX,NY);
		allocatePixels();
		
		return _fluidSolver;
	}
	
	
	//--------------------------------------------------------------
	public MSAFluidSolver setup(MSAFluidSolver f) {
		deleteFluidSolver();
		_fluidSolver = f;
		allocatePixels();
		
		return _fluidSolver;
	}
	
	//--------------------------------------------------------------
	public MSAFluidSolver getFluidSolver() {
		return _fluidSolver;
	}
	
	//--------------------------------------------------------------
	public void enableAlpha(bool b) {
		_alphaEnabled = b;
		
		if(isFluidReady()) {
			allocatePixels();
		}
	}
	
	
	//--------------------------------------------------------------
	public void allocatePixels() {
		int texWidth = _fluidSolver.getWidth()-2;
		int texHeight =_fluidSolver.getHeight()-2;
		_pixels = new Color[texWidth * texHeight];
	}
	
	
	//--------------------------------------------------------------
	public void reset() {
		if(!isFluidReady()) {
			Debug.Log("reset() - Fluid not ready\n");
			return;
		}
		_fluidSolver.reset();
	}
	
	//--------------------------------------------------------------
	public void update() {
		if(!isFluidReady()) {
			Debug.Log("updateFluid() - Fluid not ready\n");
			return;
		}
		_fluidSolver.update();
	}
	
	
	//--------------------------------------------------------------
	public void setDrawMode(DrawMode newDrawMode) {
		drawMode = newDrawMode;
		if(drawMode < 0) drawMode = (DrawMode.kDrawCount-1);
		else if(drawMode >= DrawMode.kDrawCount) drawMode = (DrawMode)0;
	}
	
	//--------------------------------------------------------------
	public void incDrawMode() {
		setDrawMode((DrawMode)((int)drawMode + 1));
	}
	
	//--------------------------------------------------------------
	public void decDrawMode() {
		setDrawMode((DrawMode)(drawMode - 1));
	}
	
	//--------------------------------------------------------------
	public DrawMode getDrawMode() {
		return drawMode;
	}
	
	//--------------------------------------------------------------
	public string getDrawModeName() {
		return getDrawModeTitles()[(int)drawMode];
	}

	//--------------------------------------------------------------
	public void draw()  {
		if(enabled == false) return;
		
		switch(drawMode) {
		case DrawMode.kDrawColor:
			drawColor();
			break;
			
		case DrawMode.kDrawMotion:
			drawMotion();
			break;
			
		case DrawMode.kDrawSpeed:
			drawSpeed();
			break;
			
		case DrawMode.kDrawVectors:
			drawVectors();
			break;
			
		default:
			break;
			
		}
	}
	
	
	//--------------------------------------------------------------
	public void drawColor(bool withAlpha = false)  {
		if(enabled == false) return;
		
		int fw = _fluidSolver.getWidth();
		int fh = _fluidSolver.getHeight();
		
		Vec2f vel;
		Color color;
		int index = 0;
		for(int j=1; j < fh-1; j++) {
			for(int i=1; i < fw-1; i++) {
				_fluidSolver.getInfoAtCell(i, j, out vel, out color);
				float r = min(color.r * brightness, 1f);
				float g = min(color.g * brightness, 1f);
				float b = min(color.b * brightness, 1f);
				if(doInvert) {
					r = 1f - r;
					g = 1f - g;
					b = 1f - b;
				}
				_pixels[index].r = r;
				_pixels[index].g = g;
				_pixels[index].b = b;
				
				if(_alphaEnabled) _pixels[index].a = withAlpha ? max(b, max(r, g)) : 1f;
				index++;
			}
		}
	}
	
	
	
	//--------------------------------------------------------------
	public void drawMotion(bool withAlpha = false)  {
		if(enabled == false) return;
		
		int fw = _fluidSolver.getWidth();
		int fh = _fluidSolver.getHeight();
		
		Vec2f vel;
		Color color;
		int index = 0;
		for(int j=1; j < fh-1; j++) {
			for(int i=1; i < fw-1; i++) {
				_fluidSolver.getInfoAtCell(i, j, out vel,out color);
				float speed2 = fabs(vel.x) * fw + fabs(vel.y) * fh;
				int speed = (int)min(speed2 * brightness, 1f);
				_pixels[index].r = min(fabs(vel.x) * fw * brightness, 1f);
				_pixels[index].g = min(fabs(vel.y) * fh * brightness, 1f);
				_pixels[index] .b = 0;
				
				if(_alphaEnabled) _pixels[index].a = withAlpha ? speed : 1f;
				index++;
			}
		}
	}
	
	
	//--------------------------------------------------------------
	public void drawSpeed(bool withAlpha = false)  {
		if(enabled == false) return;

		int fw = _fluidSolver.getWidth();
		int fh = _fluidSolver.getHeight();
		
		Vec2f vel;
		Color color;
		int index = 0;
		for(int j=1; j < fh-1; j++) {
			for(int i=1; i < fw-1; i++) {
				_fluidSolver.getInfoAtCell(i, j, out vel,out color);
				float speed2 = fabs(vel.x) * fw + fabs(vel.y) * fh;
				int speed = (int)min(speed2 * brightness, 1f);
				_pixels[index].r = speed;
				_pixels[index].g = speed;
				_pixels[index].b = speed;
				
				if(_alphaEnabled) _pixels[index].a = withAlpha ? speed : 1f;
				index++;
			}
		}
	}
	
	
	//--------------------------------------------------------------
	public void drawVectors()   {
		if(enabled == false) return;
		
		int fw = _fluidSolver.getWidth();
		int fh = _fluidSolver.getHeight();

		
		float maxVel = 5.0f/20000;
		
		Vec2f vel;
		float vt = velDrawThreshold * _fluidSolver.getInvWidth() * _fluidSolver.getInvHeight();
		vt *= vt;
		int index;
		for (int j=0; j<fh-2; j+=vectorSkipCount+1 ){
			for (int i=0; i<fw-2; i+=vectorSkipCount+1 ){
				vel = _fluidSolver.getVelocityAtCell(i+1, j+1);
				float d2 = vel.SqrMagnitude();
				if(d2>vt) {
					if(d2 > maxVel * maxVel) {
						float mult = maxVel * maxVel/ d2;
						//				float mult = (d2 - maxVel * maxVel) / d2;
						vel.x *= mult;
						vel.y *= mult;
					}
					vel *= velDrawMult * 50000;

					index = j *(fw - 1) + i;

					_pixels[i].r = vel.x;
					_pixels[i].g = vel.y;
					_pixels[i].b = 0f;
					_pixels[i].a = 1f;
				}
			}
		}
	}
	
	
	
	//--------------------------------------------------------------
	void deleteFluidSolver() {
		//	Debug.Log("deleteFluidSolver()\n");
		if(_fluidSolver != null && _didICreateTheFluid) {
			_fluidSolver = null;
			_pixels = null;
		}
	}
	
	//--------------------------------------------------------------
	bool isFluidReady() {
		if(_fluidSolver == null) {
			Debug.Log("isFluidReady() - No fluid solver\n");
			return false;
		}
		
		if(!_fluidSolver.isInited()) {
			Debug.Log("isFluidReady() - Fluid solver not initialized yet, call init()\n");
			return false;
		}
		
		return true;
	}
}

