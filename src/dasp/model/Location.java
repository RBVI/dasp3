
package dasp.model;

public class Location {
	double x;
	double y;
	double z;

	public Location (double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public double getX() { return x; }
	public double getY() { return y; }
	public double getZ() { return z; }
}
