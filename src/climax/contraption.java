/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package climax;

import java.awt.Graphics2D;

/**
 *
 * @author pchr
 */
public interface contraption {
    public int getID();
    public void SelfPortrait(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs);
    public void Motion(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs, int step, double scale);
    public double maximum_X_coordinate();
    public double minimum_X_coordinate();
    public double maximum_Y_coordinate();
    public double minimum_Y_coordinate();
}
