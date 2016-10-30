using System.Collections;
using UnityEngine;
using Vec2f = UnityEngine.Vector2;
using Vec3f = UnityEngine.Vector3;
using Color = UnityEngine.Color;

public partial class MSAFluidSolver : MSAFluidSolverUtililty
{
	// get fluid cell index for cell coordinates or normalized position
	//-------- get index
	public int getIndexForCell(int i, int j)  {
		i = clamp(i, 1, _NX);
		j = clamp(j, 1, _NY);
		return FLUID_IX(i, j);
	}
	
	public int getIndexForPos( Vec2f pos)  {
		return getIndexForCell((int)floor(pos.x * width), (int)floor(pos.y * height));
	}

	//-------- get info
	public void getInfoAtIndex(int index,out Vec2f vel, out Color color)  {
		vel = getVelocityAtIndex(index);
		color = getColorAtIndex(index);
	}
	
	public void getInfoAtCell(int i, int j,out Vec2f vel,out Color color)  {
		getInfoAtIndex(getIndexForCell(i, j),out vel,out color);
	}
	
	
	public void getInfoAtPos( Vec2f pos, out Vec2f vel, out Color color)  {
		getInfoAtIndex(getIndexForPos(pos),out vel,out color);
	}
	
	
	//-------- get velocity
	public Vec2f getVelocityAtIndex(int index)  {
		return uv[index];
	}
	
	public Vec2f getVelocityAtCell(int i, int j)  {
		return getVelocityAtIndex(getIndexForCell(i, j));
	}
	
	public Vec2f getVelocityAtPos( Vec2f pos)  {
		return getVelocityAtIndex(getIndexForPos(pos));
	}
	
	
	//-------- get color
	public Color getColorAtIndex(int index)  {
		if(doRGB) { 
			return new Color(this.color[index].x, this.color[index].y, this.color[index].z);
		} else {
			return new Color(density[index], density[index], density[index]);
		}
	}
	
	public Color getColorAtCell(int i, int j)  {
		return getColorAtIndex(getIndexForCell(i, j));
	}
	
	public Color getColorAtPos( Vec2f pos)  {
		return getColorAtIndex(getIndexForPos(pos));
	}
	
	
	//-------- add force
	public void addForceAtIndex(int index,  Vec2f force) {
		uv[index] += force;
	}
	
	public void addForceAtCell(int i, int j,  Vec2f force) {
		addForceAtIndex(getIndexForCell(i, j), force);
	}
	
	public void addForceAtPos( Vec2f pos,  Vec2f force) {
		addForceAtIndex(getIndexForPos(pos), force);
	}
	
	
	//-------- add color
	public void addColorAtIndex(int index,  Vec3f color) {
		if(doRGB) {
			colorOld[index] += new Vec3f(color.x, color.y, color.z);
		} else {
			density[index] += color.x;
		}
	}
	
	public void addColorAtCell(int i, int j,  Vec3f color) {
		addColorAtIndex(getIndexForCell(i, j), color);
	}
	
	public void addColorAtPos( Vec2f pos,  Vec3f color) {
		addColorAtIndex(getIndexForPos(pos), color);
	}

	public void addSource (Vec2f[] x, Vec2f[] x0) {
		for(int i = _numCells-1; i >=0; --i) {
			x[i] += x0[i] * deltaT;
		}
	}

	public void addSource (Vec3f[] x, Vec3f[] x0) {
		for(int i = _numCells-1; i >=0; --i) {
			x[i] += x0[i] * deltaT;
		}
	}

	public void addSource (float[] x, float[] x0) {
		for(int i = _numCells-1; i >=0; --i) {
			x[i] += x0[i] * deltaT;
		}
	}

	public int getIndexForNormalizedPosition(float x, float y) {
		return getIndexForCell((int)Mathf.Floor(x * (_NX)), (int)Mathf.Floor(y * (_NY)));
	}
}
