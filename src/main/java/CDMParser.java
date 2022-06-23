import com.google.gson.Gson;
import org.hipparchus.analysis.function.Abs;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.Map;
import java.util.Scanner;

public class CDMParser {
    public static void main(String[] arg) throws FileNotFoundException {

        DataLoader loader = new DataLoader(); //creates an instance of the loader type
        loader.load(); // loads the orekit-data file to get constant parameters like tai-utc, etc.

        //define the reference frames
        Frame eme2000 = FramesFactory.getEME2000();
        Frame gcrf = FramesFactory.getGCRF();
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010,false);

        String spPath = "/Users/connergrey/Documents/SP VECTORS/vectors_22115/scratch/SP_VEC";
        SPVecReader spVecReader = new SPVecReader();

        StringBuffer buf = new StringBuffer();
        TimeScale utc = TimeScalesFactory.getUTC();
        PrintWriter covWriter = new PrintWriter("covPropTesting.csv");

        buf.append("Creation Date");
        buf.append(",");

        buf.append("Last OD Date");
        buf.append(",");

        buf.append("Norad ID");
        buf.append(",");

        buf.append("TCA");
        buf.append(",");

        buf.append("Days from sp epoch");
        buf.append(",");

        buf.append("Days from last od");
        buf.append(",");

        buf.append("X");
        buf.append(",");
        buf.append("Y");
        buf.append(",");
        buf.append("Z");
        buf.append(",");
        buf.append("VX");
        buf.append(",");
        buf.append("VY");
        buf.append(",");
        buf.append("VZ");
        buf.append(",");

        buf.append("RR");
        buf.append(",");
        buf.append("IR");
        buf.append(",");
        buf.append("II");
        buf.append(",");
        buf.append("CR");
        buf.append(",");
        buf.append("CI");
        buf.append(",");
        buf.append("CC");
        buf.append(",");
        buf.append("CRDOT_R");
        buf.append(",");
        buf.append("CRDOT_T");
        buf.append(",");
        buf.append("CRDOT_N");
        buf.append(",");
        buf.append("CRDOT_RDOT");
        buf.append(",");
        buf.append("CTDOT_R");
        buf.append(",");
        buf.append("CTDOT_T");
        buf.append(",");
        buf.append("CTDOT_N");
        buf.append(",");
        buf.append("CTDOT_RDOT");
        buf.append(",");
        buf.append("CTDOT_TDOT");
        buf.append(",");
        buf.append("CNDOT_R");
        buf.append(",");
        buf.append("CNDOT_T");
        buf.append(",");
        buf.append("CNDOT_N");
        buf.append(",");
        buf.append("CNDOT_RDOT");
        buf.append(",");
        buf.append("CNDOT_TDOT");
        buf.append(",");
        buf.append("CNDOT_NDOT");


        buf.append("\n");

        int i = 0;
        int j = 0;
        File cdmFolder = new File("/Users/connergrey/cdms");
        for (File cdmSubFolder: cdmFolder.listFiles()) {
            if (cdmSubFolder.toString().compareTo(cdmFolder.toString() + "/.DS_Store") == 0){
                continue;
            }

            for (File cdm : cdmSubFolder.listFiles()) {
                if (cdm.toString().compareTo(cdmSubFolder.toString() + "/.DS_Store") == 0){
                    continue;
                }

                // addd loop through all json files
                Scanner scan = new Scanner(cdm);
                //Scanner scan = new Scanner(new File("275029489.json"));

                String[] pts = scan.nextLine().split("\\{");
                String[] pts2 = pts[2].split("\\}");
                String jsonStr = "{" + pts2[0] + "}";
                //System.out.println(jsonStr);

                PrintWriter pw = new PrintWriter("test.json");
                pw.write(jsonStr);
                pw.close();

                try {
                    // create Gson instance
                    Gson gson = new Gson();
                    // create a reader
                    Reader reader = Files.newBufferedReader(Paths.get("test.json"));
                    // convert JSON file to map
                    Map<String, String> map = gson.fromJson(reader, Map.class);

                    AbsoluteDate creation = new AbsoluteDate(map.get("CREATION_DATE"), utc);
                    AbsoluteDate tca = new AbsoluteDate(map.get("TCA"), utc);


                    if (map.get("SAT1_EPHEMERIS_NAME").equals("NONE")) {

                        if(Integer.parseInt( map.get("SAT1_OBJECT_DESIGNATOR") ) > 52260 ){
                            // print map entries
//                            for (Map.Entry<?, ?> entry : map.entrySet()) {
//                                System.out.println(entry.getKey() + "=" + entry.getValue());
//                            }
                            continue;
                        }

                        AbsoluteDate lastOD = new AbsoluteDate(map.get("SAT1_TIME_LASTOB_END"), utc);
                        double durationFromLastOD = tca.durationFrom(lastOD) / 86400;
                        double creationFromOD = creation.durationFrom(lastOD) / 86400;

                        //read the sp vector
                        spVecReader.readSPVec(spPath,Integer.parseInt(map.get("SAT1_OBJECT_DESIGNATOR")));
                        AbsoluteDate spEpoch = spVecReader.getDate();
                        double spEpochtoTCA = tca.durationFrom(spEpoch) / 86400;
                        double spEpochtoLastOD = lastOD.durationFrom(spEpoch)/ 86400;

                        //read the frame of reference
                        String refFrame = map.get("SAT1_REF_FRAME");
                        //define the frame transformations
                        Transform trans = null;
                        switch (refFrame) {
                            case "ITRF" -> {
                                trans = itrf.getTransformTo(eme2000, tca);
                            }
                            case "GCRF" -> {
                                trans = gcrf.getTransformTo(eme2000, tca);
                            }
                            case "EME2000" -> {
                                trans = eme2000.getTransformTo(eme2000, tca);
                            }
                            default -> System.out.println("error");
                        }

                        // read the position and velocity into vectors

                        PVCoordinates pvFrame = new PVCoordinates(
                                new Vector3D( Double.parseDouble(map.get("SAT1_X")),Double.parseDouble(map.get("SAT1_Y")),Double.parseDouble(map.get("SAT1_Z")) ),
                                new Vector3D( Double.parseDouble(map.get("SAT1_X_DOT")),Double.parseDouble(map.get("SAT1_Y_DOT")),Double.parseDouble(map.get("SAT1_Z_DOT")) ));

                        PVCoordinates pvEME = trans.transformPVCoordinates(pvFrame);

                        buf.append(creation.toString());
                        buf.append(",");

                        buf.append(lastOD.toString());
                        buf.append(",");

                        buf.append(map.get("SAT1_OBJECT_DESIGNATOR"));
                        buf.append(",");

                        buf.append(tca.toString());
                        buf.append(",");

                        buf.append(String.valueOf(spEpochtoTCA));
                        buf.append(",");

                        buf.append(String.valueOf(durationFromLastOD));
                        buf.append(",");

                        buf.append(pvEME.getPosition().getX());
                        buf.append(",");
                        buf.append(pvEME.getPosition().getY());
                        buf.append(",");
                        buf.append(pvEME.getPosition().getZ());

                        buf.append(",");
                        buf.append(pvEME.getVelocity().getX());
                        buf.append(",");
                        buf.append(pvEME.getVelocity().getY());
                        buf.append(",");
                        buf.append(pvEME.getVelocity().getZ());
                        buf.append(",");

                        buf.append(map.get("SAT1_CR_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CN_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CN_T"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CN_N"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CRDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CRDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CRDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CRDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CTDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CTDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CTDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CTDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CTDOT_TDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_TDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT1_CNDOT_NDOT"));
                        buf.append(",");

                        buf.append("\n");
                        covWriter.write(buf.toString());
                        buf.delete(0, buf.length());

                        i++;
                    }else{
                        j++;
                    }

                    if (map.get("SAT2_EPHEMERIS_NAME").equals("NONE")) {
                        //int sat1No = Integer.parseInt(map.get("SAT1_OBJECT_DESIGNATOR"));

                        // removes the ones that have zeros for the covariance matrix and report MD only
                        // all seem to have null covariance and null LASTOB END TIME
                        // throws errors if I dont remove
                        // all seem to be in the 80000+ norad id range
                        if(Integer.parseInt( map.get("SAT2_OBJECT_DESIGNATOR") ) > 52260 ){
                            // print map entries
//                            for (Map.Entry<?, ?> entry : map.entrySet()) {
//                                System.out.println(entry.getKey() + "=" + entry.getValue());
//                            }
                            continue;
                        }
                        AbsoluteDate lastOD = new AbsoluteDate(map.get("SAT2_TIME_LASTOB_END"), utc);
                        double durationFromLastOD = tca.durationFrom(lastOD) / 86400;
                        double creationFromOD = creation.durationFrom(lastOD) / 86400;

                        //read the sp vector
                        spVecReader.readSPVec(spPath,Integer.parseInt(map.get("SAT2_OBJECT_DESIGNATOR")));
                        AbsoluteDate spEpoch = spVecReader.getDate();
                        double spEpochtoTCA = tca.durationFrom(spEpoch) / 86400;
                        double spEpochtoLastOD = lastOD.durationFrom(spEpoch) / 86400;

                        //read the frame of reference
                        String refFrame = map.get("SAT2_REF_FRAME");
                        //define the frame transformations
                        Transform trans = null;
                        switch (refFrame) {
                            case "ITRF" -> {
                                trans = itrf.getTransformTo(eme2000, tca);
                            }
                            case "GCRF" -> {
                                trans = gcrf.getTransformTo(eme2000, tca);
                            }
                            case "EME2000" -> {
                                trans = eme2000.getTransformTo(eme2000, tca);
                            }
                            default -> System.out.println("error");
                        }

                        // read the position and velocity into vectors

                        PVCoordinates pvFrame = new PVCoordinates(
                                new Vector3D( Double.parseDouble(map.get("SAT2_X")),Double.parseDouble(map.get("SAT2_Y")),Double.parseDouble(map.get("SAT2_Z")) ),
                                new Vector3D( Double.parseDouble(map.get("SAT2_X_DOT")),Double.parseDouble(map.get("SAT2_Y_DOT")),Double.parseDouble(map.get("SAT2_Z_DOT")) ));

                        PVCoordinates pvEME = trans.transformPVCoordinates(pvFrame);

                        buf.append(creation.toString());
                        buf.append(",");

                        buf.append(lastOD.toString());
                        buf.append(",");

                        buf.append(map.get("SAT2_OBJECT_DESIGNATOR"));
                        buf.append(",");

                        buf.append(tca.toString());
                        buf.append(",");

                        buf.append(String.valueOf(spEpochtoTCA));
                        buf.append(",");

                        buf.append(String.valueOf(durationFromLastOD));
                        buf.append(",");

                        buf.append(pvEME.getPosition().getX());
                        buf.append(",");
                        buf.append(pvEME.getPosition().getY());
                        buf.append(",");
                        buf.append(pvEME.getPosition().getZ());
                        buf.append(",");
                        buf.append(pvEME.getVelocity().getX());
                        buf.append(",");
                        buf.append(pvEME.getVelocity().getY());
                        buf.append(",");
                        buf.append(pvEME.getVelocity().getZ());
                        buf.append(",");

                        buf.append(map.get("SAT2_CR_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CN_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CN_T"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CN_N"));
                        buf.append(",");

                        buf.append(map.get("SAT2_CRDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CRDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CRDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CRDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CTDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CTDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CTDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CTDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CTDOT_TDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_R"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_T"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_N"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_RDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_TDOT"));
                        buf.append(",");
                        buf.append(map.get("SAT2_CNDOT_NDOT"));

                        buf.append("\n");
                        covWriter.write(buf.toString());
                        buf.delete(0, buf.length());
                        i++;
                    }else{
                        j++;
                    }

                    // close reader
                    reader.close();

                } catch (Exception ex) {
                    ex.printStackTrace();
                }


            /*if(i>10)
                break;*/
            }
        }
        covWriter.close();
        System.out.println("Included " + i + " entries");
        System.out.println("Excluded " + j + " entries");


    }
}
