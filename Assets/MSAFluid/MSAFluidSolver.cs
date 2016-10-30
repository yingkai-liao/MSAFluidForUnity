using UnityEngine;
using System.Collections;
using Vec2f = UnityEngine.Vector2;
using Vec3f = UnityEngine.Vector3;

public partial class MSAFluidSolver : MSAFluidSolverUtililty
{
	// do not change these values, you can override them using the solver methods
	public const	int FLUID_DEFAULT_NX					= 100;
	public const	int FLUID_DEFAULT_NY					= 100;
	public const	float FLUID_DEFAULT_DT					= 0.04f;	//Maa	25fps
	public const	float FLUID_DEFAULT_VISC				= 0.0001f;
	public const	float FLUID_DEFAULT_COLOR_DIFFUSION		= 0.0f;
	public const	float FLUID_DEFAULT_FADESPEED				= 0.03f;
	public const	int FLUID_DEFAULT_SOLVER_ITERATIONS	= 10;

	public int	FLUID_IX(int i,int j){return (i) + (_NX + 2)  *(j);}

		// allocate an array large enough to hold information for u, v, r, g, OR b
	public float[] alloc()	{ return new float[_numCells];	}
		
		
	public float[] density, densityOld;		// used if not RGB
	public Vec3f[]	color, colorOld;			// used for RGB
	public Vec2f[]	uv, uvOld;
		
	public float[]	curl;
		
	public bool	doRGB;				// for monochrome, update only density
	public bool	doVorticityConfinement;
	public int		solverIterations;
		
	public float	colorDiffusion;
	public float	viscocity;
	public float	fadeSpeed;
	public float	deltaT;
	public bool	wrap_x;
	public bool	wrap_y;
		
	//protected:
			
	protected 	float width;
	protected float height;
	protected float invWidth;
	protected float invHeight;
		
	protected int		_NX, _NY, _numCells;
	protected float	_invNX, _invNY, _invNumCells;
	protected bool	_isInited;
	protected 	float[] _tmp;
		
	protected 	float	_avgDensity;			// this will hold the average color of the last frame (how full it is)
	protected 	float	_uniformity;			// this will hold the _uniformity of the last frame (how uniform the color is);
	protected 	float	_avgSpeed;
}
			
			