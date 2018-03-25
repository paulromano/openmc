from copy import copy

import openmc


class RightCircularCylinder(openmc.Surface):
    def __init__(self, center_base, height, radius, axis='z',
                 boundary_type='transmission'):
        kwargs = {'boundary_type': boundary_type}
        cx, cy, cz = center_base
        if axis == 'x':
            self.cyl = openmc.XCylinder(y0=cy, z0=cz, R=radius, **kwargs)
            self.bottom = openmc.XPlane(x0=cx, **kwargs)
            self.top = openmc.XPlane(x0=cx + height, **kwargs)
        elif axis == 'y':
            self.cyl = openmc.YCylinder(x0=cx, z0=cz, R=radius, **kwargs)
            self.bottom = openmc.YPlane(y0=cy, **kwargs)
            self.top = openmc.YPlane(y0=cy + height, **kwargs)
        elif axis == 'z':
            self.cyl = openmc.ZCylinder(x0=cx, y0=cy, R=radius, **kwargs)
            self.bottom = openmc.ZPlane(z0=cz, **kwargs)
            self.top = openmc.ZPlane(z0=cz + height, **kwargs)

    def __neg__(self):
        return -self.cyl & +self.bottom & -self.top

    def __pos__(self):
        return +self.cyl | -self.bottom | +self.top

    def evaluate(self, point):
        raise NotImplementError('Macrobodies do not have a surface equation.')

    def translate(self, vector):
        surf = copy(self)
        surf.cyl = surf.cyl.translate(vector)
        surf.bottom = surf.bottom.translate(vector)
        surf.top = surf.top.translate(vector)
        return surf


class Box(openmc.Surface):
    def __init__(self, v, a1, a2, a3, boundary_type='transmission'):
        kwargs = {'boundary_type': boundary_type}
        vx, vy, vz = v
        a1x, a1y, a1z = a1
        a2x, a2y, a2z = a2
        a3x, a3y, a3z = a3

        # Only support boxes with axis-aligned vectors
        if any(x != 0.0 for x in (a1y, a1z, a2x, a2z, a3x, a3y)):
            raise NotImplementedError('Box macrobody with non-axis-aligned '
                                      'vector not supported.')

        # Determine each side of the box
        if a1x > 0:
            xmin, xmax = vx, vx + a1x
        else:
            xmin, xmax = vx + a1x, vx
        if a2y > 0:
            ymin, ymax = vy, vy + a2y
        else:
            ymin, ymax = vy + a2y, vy
        if a3z > 0:
            zmin, zmax = vz, vz + a3z
        else:
            zmin, zmax = vz + a3z, vz

        # Create surfaces
        self.xmin = openmc.XPlane(x0=xmin, **kwargs)
        self.xmax = openmc.XPlane(x0=xmax, **kwargs)
        self.ymin = openmc.YPlane(y0=ymin, **kwargs)
        self.ymax = openmc.YPlane(y0=ymax, **kwargs)
        self.zmin = openmc.ZPlane(z0=zmin, **kwargs)
        self.zmax = openmc.ZPlane(z0=zmax, **kwargs)

    def __neg__(self):
        return (+self.xmin & -self.xmax &
                +self.ymin & -self.ymax &
                +self.zmin & -self.zmax)

    def __pos__(self):
        return (-self.xmin | +self.xmax |
                -self.ymin | +self.ymax |
                -self.zmin | +self.zmax)

    def evaluate(self, point):
        raise NotImplementError('Macrobodies do not have a surface equation.')

    def translate(self, vector):
        surf = copy(self)
        surf.xmin = surf.xmin.translate(vector)
        surf.xmax = surf.xmax.translate(vector)
        surf.ymin = surf.ymin.translate(vector)
        surf.ymax = surf.ymax.translate(vector)
        surf.zmin = surf.zmin.translate(vector)
        surf.zmax = surf.zmax.translate(vector)
        return surf
