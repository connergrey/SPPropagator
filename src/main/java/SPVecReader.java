import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.TimeStampedPVCoordinates;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Calendar;
import java.util.Scanner;

public class SPVecReader {

    private int degOrd;
    private double cD;
    private double cR;
    private AbsoluteDate date;
    private int satNo;
    private TimeStampedPVCoordinates tspvTEME;

    public void readSPVec(String path, int satNo){

        // read the SP vector that the conjunction was derived from
        String satNoString;
        if(Integer.toString(satNo).length()==1){
            satNoString = "0000"+Integer.toString(satNo);
        }else if(Integer.toString(satNo).length()==2){
            satNoString = "000"+Integer.toString(satNo);
        }else if(Integer.toString(satNo).length()==3){
            satNoString = "00"+Integer.toString(satNo);
        }else if(Integer.toString(satNo).length()==4){
            satNoString = "0"+Integer.toString(satNo);
        }else{
            satNoString = Integer.toString(satNo);
        }
        //pull the sp vector
        String filePath = path + "/" + satNoString.substring(0,2) + "/" + satNoString + ".txt";

        readSPVec(new File(filePath));
    }

    public void readSPVec(String stringFileName){
        readSPVec(new File(stringFileName));
    }

    public void readSPVec(File fName){

        try{
            //System.out.println(fName);
            Scanner scan = new Scanner(fName);

            // GET SAT NO
            scan.next();
            scan.next();
            satNo = scan.nextInt();

            // GET DATE
            scan.nextLine();
            scan.next();
            scan.next();
            scan.next();
            int year = scan.nextInt();
            int dayOfYear = scan.nextInt();

            scan.next();
            scan.next();
            String time = scan.next();
            int hour = Integer.parseInt(time.substring(0,2));
            int min = Integer.parseInt(time.substring(3,5));

            double sec = 0;
            if(time.length() == 6) {
                sec = scan.nextDouble();
            }else if(time.length() == 12){
                sec = Double.parseDouble(time.substring(6,time.length()));
            }else{
                System.out.println("theres still an error");
            }

            Calendar calendar = Calendar.getInstance();
            calendar.set(Calendar.YEAR, year);
            calendar.set(Calendar.DAY_OF_YEAR, dayOfYear);
            int month = calendar.get(Calendar.MONTH) + 1; // returns 0 - 11
            int day = calendar.get(Calendar.DAY_OF_MONTH);

            TimeScale utc = TimeScalesFactory.getUTC();
            date = new AbsoluteDate(year,month,day,hour,min,sec,utc);

            scan.nextLine();
            scan.next();
            scan.next();
            scan.next();
            Vector3D pos = new Vector3D(scan.nextDouble(),scan.nextDouble(),scan.nextDouble()).scalarMultiply(1000);
            scan.nextLine();
            scan.next();
            scan.next();
            scan.next();
            Vector3D vel = new Vector3D(scan.nextDouble(),scan.nextDouble(),scan.nextDouble()).scalarMultiply(1000);
            scan.nextLine();
            tspvTEME = new TimeStampedPVCoordinates(date,pos,vel);

            // GET GEOPOTENTIAL
            scan.next();
            scan.next();
            String[] par = scan.next().split("Z,");
            degOrd = Integer.parseInt( par[0] );

            // GET DRAG
            scan.next();
            String dragModel = scan.next();

            scan.nextLine();
            scan.nextLine();

            // GET COEFFS
            scan.next();
            scan.next();
            scan.next();
            cD = Double.parseDouble( scan.next() );

            scan.nextLine();
            scan.next();
            scan.next();
            scan.next();
            scan.next();
            scan.next();
            cR = Double.parseDouble( scan.next() );

            scan.close();


        }catch(FileNotFoundException e){
            e.printStackTrace();
        }
    }

    public TimeStampedPVCoordinates getInitialState(Frame frame){

        Frame teme = FramesFactory.getTEME();
        Transform temeToFrame = teme.getTransformTo(frame,date);
        return temeToFrame.transformPVCoordinates(tspvTEME);

    }



    public double getcD() {
        return cD;
    }

    public double getcR() {
        return cR;
    }

    public int getDegOrd() {
        return degOrd;
    }

    public AbsoluteDate getDate() {
        return date;
    }

    public int getSatNo() {
        return satNo;
    }
}
