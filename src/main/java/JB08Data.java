import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class JB08Data implements JB08InputParams{

    List<AbsoluteDate> dataDates = new ArrayList<>();
    int n = 200;
    double[] f10 = new double[n];
    double[] f54 = new double[n];
    double[] s10 = new double[n];
    double[] s54 = new double[n];
    double[] m10 = new double[n];
    double[] m54 = new double[n];
    double[] y10 = new double[n];
    double[] y54 = new double[n];
    double[] dtc = new double[n];


    public void readData(String fileName) {
        TimeScale utc = TimeScalesFactory.getUTC();
        // Read in 3hr interval data
        File file = new File(fileName);

        try {
            Scanner scan = new Scanner(file);
            scan.nextLine();
            int i = 0;
            while (scan.hasNextLine()) {
                String str = scan.nextLine();
                String[] splitStr = str.split(",");
                dataDates.add(new AbsoluteDate(splitStr[0], utc));
                f10[i] = Double.parseDouble(splitStr[1]);
                f54[i] = Double.parseDouble(splitStr[2]);
                s10[i] = Double.parseDouble(splitStr[3]);
                s54[i] = Double.parseDouble(splitStr[4]);
                m10[i] = Double.parseDouble(splitStr[5]);
                m54[i] = Double.parseDouble(splitStr[6]);
                y10[i] = Double.parseDouble(splitStr[7]);
                y54[i] = Double.parseDouble(splitStr[8]);
                dtc[i] = Double.parseDouble(splitStr[9]);
                i++;
            }
        } catch (
                FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Override
    public AbsoluteDate getMinDate() {
        return null;
    }

    @Override
    public AbsoluteDate getMaxDate() {
        return null;
    }

    @Override
    public double[] getData(AbsoluteDate check) {

        // Find index for date
        int ind = -1;
        for (int x = 0; x < dataDates.size(); x++) {
            if (check.compareTo(dataDates.get(x)) >= 0 && check.compareTo(dataDates.get(x + 1)) < 0) {
                ind = x;
                break;
            }
        }
        if(ind == -1){
            System.out.println("Date not found");
        }


        return new double[]{ f10[ind] , f54[ind] , s10[ind] , s54[ind] , m10[ind] , m54[ind] ,
                y10[ind] , y54[ind] , dtc[ind]
        };

    }
}
