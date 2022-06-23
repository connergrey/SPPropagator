import org.orekit.bodies.BodyShape;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.models.earth.atmosphere.NRLMSISE00;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;

import java.io.File;
import java.util.Locale;

public class JB08Example {

    public void main(String[] arg){

        //loads constant data file info
        final File home = new File(System.getProperty("user.home"));
        final File orekitData = new File(home, "orekit-data");
        if (!orekitData.exists()) {
            System.err.format(Locale.US, "Failed to find %s folder%n", orekitData.getAbsolutePath());
            System.err.format(Locale.US, "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
                    "orekit-data-master.zip", "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
                    home.getAbsolutePath());
            System.exit(1);
        }
        DataContext.
                getDefault().
                getDataProvidersManager().
                addProvider(new DirectoryCrawler(orekitData));
        // end configure

        String solarInd = "SOLFSMYDST_2022056.csv"; // solar weather data file. Is a formatted version of a UDL output command
        //comes from Space Environment Technologies, from the SGI query I believe...

        //create instance of jb08 data
        JB08Data jb08data = new JB08Data();
        jb08data.readData(solarInd); // reads from one file formed from UDL data

        //create instance of jb08 atmos model
        CelestialBody sun = CelestialBodyFactory.getSun();
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010,false);
        BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, itrf);
        TimeScale utc = TimeScalesFactory.getUTC();
        Atmosphere jb08 = new JB08(jb08data, sun, earth, utc);

        DragSensitive dragSens = new IsotropicDrag(1, 0.1);
        ForceModel drag = new DragForce(jb08, dragSens);

    }
}
