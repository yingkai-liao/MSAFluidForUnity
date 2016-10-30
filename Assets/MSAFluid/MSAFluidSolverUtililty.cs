using UnityEngine;
using System.Collections;

public class Rand
{
	public static float randFloat()
	{
		float r = Random.Range(0f,1f);
		return r;
	}
}

public class MSAFluidSolverUtililty 
{
	public int clamp(int val, int min,int max)
	{
		return Mathf.Clamp(val,min,max);
	}

	public float floor(float val)
	{
		return Mathf.FloorToInt(val);
	}

	public float sqrt(float val)
	{
		return Mathf.Sqrt(val);
	}

	public float fabs(float val)
	{
		return Mathf.Abs(val);
	}

	public float min(float v1,float v2)
	{
		return Mathf.Min(v1,v2);
	}

	public float max(float v1,float v2)
	{
		return Mathf.Max(v1,v2);
	}

	public void SWAP<T>(ref T a,ref T b)						
	{
		T tmp = b; b = a; a = tmp; 
	}
}