
import general.utils.DataUtils;
import general.utils.RDataUtils;
import org.rosuda.REngine.Rserve.RConnection;

public class UploadBean implements Serializable {

    private final ApplicationBean1 ab = (ApplicationBean1) DataUtils.findBean("applicationBean1");
    private final SessionBean1 sb = (SessionBean1) DataUtils.findBean("sessionBean1");
	
	private UploadedFile file;

   
    public String handleomicsTextFileUpload() {
        
        if (file == null || file.getSize() <= 0 || samplefilet == null || samplefilet.getSize() <= 0) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Please check your omics data!"));
            return null;
        } else{
            sb.doLogin("o2pls"); //analType
            sb.setModuleType("omics");
            sb.setDataFormat("text");
            
            RConnection RC = sb.getRConnection();
            String homeDir = sb.getCurrentUser().getHomeDir();
            
            String fileName = file.getFileName(); 
            if (!DataUtils.uploadFile(file, homeDir)) {
                return null;
            }

            String samfile = samplefilet.getFileName();
            if (!DataUtils.uploadFile(samplefilet, homeDir)) {
                return null;
            }


            int res3 = RDataUtils.readOmicsTable(RC, fileName, samfile);
            
            if (res3 == 0){
                String errMsg = RDataUtils.getCurrentMsg(RC);
                FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
                return null;                
            }else {
                String format = sb.getDataFormat();
                String datatype = sb.getModuleType();
        
                int res = RDataUtils.missingCheckData(RC, datatype);
                if (res == 1) {
                    FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_INFO, "Only Omics 1 has missing values.", "OK"));
                    return "Missing value1";
                } else if (res == 2){
                    FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_INFO, "Only Omics 2 has missing values.", "OK"));
                    return "Missing value2";
                } else if (res == 3){
                    FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_INFO, "Omics 1 and 2 have missing values.", "OK"));
                    return "Missing value12";
                } else {
                    FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_INFO, "Both data has no missing value.", "OK"));
                    if (setup16SUploadInfo(fileName, samfile) == 1) { //this part is going to the next page
                        return "Data Check";
                    }
                }
 
            }
            
        }
        return null;
    }
    
    //After missing value imputed...
    public String handleomicsTextFileUpload2() {
        RConnection RC = sb.getRConnection();
        String format = sb.getDataFormat();
        String datatype = sb.getModuleType();
            
        if (performomicsdataprocessing3() == 1) { //this part is going to the next page
            return "Data Check";
        }
        return null;
    }
    
    //After missing value imputed...
    public int performomicsdataprocessing3() {
        
        
        RConnection RC = sb.getRConnection();
        String datatype = sb.getModuleType();
        int[] res = RDataUtils.sanityCheckData3(RC, datatype);
        if (res[0] == 1) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data has been uploaded sucessfully. All required files are present for analysis.", "OK"));
            String msg = "<table>";
            if (datatype.equals("omics")) {

                msg = msg + "<tr><td align=\"left\" width=\"300\"> <b>Data type:</b></td><td align=\"left\">" + datatype + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>File format: </b></td><td align=\"left\">" + sb.getDataFormat() + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics1:</b></td><td align=\"left\"> " + res[1] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics2:</b></td><td align=\"left\"> " + res[2] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics1: </b></td><td align=\"left\"> " + res[3] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics2: </b></td><td align=\"left\"> " + res[4] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics1: </b></td><td align=\"left\"> " + res[5] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics2: </b></td><td align=\"left\"> " + res[6] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics1: </b></td><td align=\"left\">" + res[7] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics2: </b></td><td align=\"left\">" + res[8] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics1: </b></td><td align=\"left\">" + res[9] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics2: </b></td><td align=\"left\">" + res[10] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics1: </b></td><td align=\"left\">" + res[11] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics2: </b></td><td align=\"left\">" + res[12] + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample names in metadata and omics matched: </b></td><td align=\"left\">" + ((res[13]==1)?"Yes":"<font style='color:red'>No</font>") + "</td></tr>";
            } else {
                String errMsg = RDataUtils.getCurrentMsg(RC);
                FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            
            }
            msg = msg + "</table>";

            sb.setProcMsg(msg);

            return (1);
        } else {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            return (0);
        }
    }
    
    
    
    private int setup16SUploadInfo(String fileName, String samfile) {
        RConnection RC = sb.getRConnection();
        String procMsg = RDataUtils.getCurrentMsg(RC);
        FacesContext.getCurrentInstance().addMessage(null,
                new FacesMessage(FacesMessage.SEVERITY_INFO, "Upload - OK", fileName + " and " + samfile + " are uploaded and parsed out."
                        + procMsg));
        
        int res = performomicsdataprocessing();
        return (res);
    }
    
            
    public int performomicsdataprocessing() {
        
        
        RConnection RC = sb.getRConnection();
        String datatype = sb.getModuleType();
        int[] res = RDataUtils.sanityCheckData(RC, datatype);
        if (res[0] == 1) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data has been uploaded sucessfully. All required files are present for analysis.", "OK"));
            String msg = "<table>";
            if (datatype.equals("omics")) {

                msg = msg + "<tr><td align=\"left\" width=\"300\"> <b>Data type:</b></td><td align=\"left\">" + datatype + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>File format: </b></td><td align=\"left\">" + sb.getDataFormat() + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics1:</b></td><td align=\"left\"> " + res[1] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics2:</b></td><td align=\"left\"> " + res[2] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics1: </b></td><td align=\"left\"> " + res[3] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics2: </b></td><td align=\"left\"> " + res[4] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics1: </b></td><td align=\"left\"> " + res[5] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics2: </b></td><td align=\"left\"> " + res[6] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics1: </b></td><td align=\"left\">" + res[7] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics2: </b></td><td align=\"left\">" + res[8] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics1: </b></td><td align=\"left\">" + res[9] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics2: </b></td><td align=\"left\">" + res[10] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics1: </b></td><td align=\"left\">" + res[11] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics2: </b></td><td align=\"left\">" + res[12] + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample names in metadata and omics matched: </b></td><td align=\"left\">" + ((res[13]==1)?"Yes":"<font style='color:red'>No</font>") + "</td></tr>";
            } else {
                String errMsg = RDataUtils.getCurrentMsg(RC);
                FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            
            }
            msg = msg + "</table>";

            sb.setProcMsg(msg);
            return (1);
        } else {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            return (0);
        }
    }
    
    
    
    //This is when the data is already normalized and filtered well -> For example datasets 
    private int setup16SUploadInfo2(String fileName, String samfile) {
        RConnection RC = sb.getRConnection();
        String procMsg = RDataUtils.getCurrentMsg(RC);
        FacesContext.getCurrentInstance().addMessage(null,
                new FacesMessage(FacesMessage.SEVERITY_INFO, "Upload - OK", fileName + " and " + samfile + " are uploaded and parsed out."
                        + procMsg));
        
        int res = performomicsdataprocessing2();
        return (res);
    }
    
   //This is when the data is already normalized and filtered well -> For example datasets         
    public int performomicsdataprocessing2() {
        
        
        RConnection RC = sb.getRConnection();
        String datatype = sb.getModuleType();
        int[] res = RDataUtils.CheckData(RC, datatype);
        if (res[0] == 1) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data has been uploaded sucessfully. All required files are present for analysis.", "OK"));
            String msg = "<table>";
            if (datatype.equals("omics")) {

                msg = msg + "<tr><td align=\"left\" width=\"300\"> <b>Data type:</b></td><td align=\"left\">" + datatype + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>File format: </b></td><td align=\"left\">" + sb.getDataFormat() + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics1:</b></td><td align=\"left\"> " + res[1] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Sample number for omics2:</b></td><td align=\"left\"> " + res[2] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics1: </b></td><td align=\"left\"> " + res[3] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Number of features for omics2: </b></td><td align=\"left\"> " + res[4] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics1: </b></td><td align=\"left\"> " + res[5] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Total values for omics2: </b></td><td align=\"left\"> " + res[6] + "</td></tr>";
                
                
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics1: </b></td><td align=\"left\">" + res[7] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Average values per sample for omics2: </b></td><td align=\"left\">" + res[8] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics1: </b></td><td align=\"left\">" + res[9] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Minimum values per sample for omics2: </b></td><td align=\"left\">" + res[10] + "</td></tr>";
                
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics1: </b></td><td align=\"left\">" + res[11] + "</td></tr>";
                msg = msg + "<tr><td align=\"left\"> <b>Maximum values per sample for omics2: </b></td><td align=\"left\">" + res[12] + "</td></tr>";

                msg = msg + "<tr><td align=\"left\"> <b>Sample names in metadata and omics matched: </b></td><td align=\"left\">" + ((res[13]==1)?"Yes":"<font style='color:red'>No</font>") + "</td></tr>";
            } else {
                String errMsg = RDataUtils.getCurrentMsg(RC);
                FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            
            }
            msg = msg + "</table>";

            sb.setProcMsg(msg);

            return (1);
        } else {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            return (0);
        }
    }
    
    
    
    
    //Example datasets selections
    private String selectedTestData = "omics1";

    public String getSelectedTestData() {
        return selectedTestData;
    }

    public void setSelectedTestData(String selectedTestData) {
        this.selectedTestData = selectedTestData;
    }
    
    //It used to be doDefault16SAnalysis
    public String doDefaultomicsAnalysis() {
        sb.doLogin("o2pls"); //analType
 
        if (selectedTestData.equals("omics2")) {
            doDefaultDILGOMAnalysis();
        }else if (selectedTestData.equals("omics3")){
            doDefaultCANCERAnalysis();
        }else if (selectedTestData.equals("omics4")){
            doDefaultNUTRIMICEAnalysis();
        }else if (selectedTestData.equals("omics5")){
            doDefaultIBDHUMANAnalysis();
        }else { //when omics1
            doDefaultNCI60Analysis();
        }
        return "Data Check";
    }
    
    //It used to be doDefaultMammalianAnalysis
    private void doDefaultNCI60Analysis() {
        sb.setModuleType("omics");
        sb.setDataFormat("text");
        

        String inpath = ab.getTestmrnaDataPath();
        String inpatht = ab.getTestproteinDataPath();
        
        
        String omics1 = DataUtils.getJustFileName(inpath);
        String omics2 = DataUtils.getJustFileName(inpatht);
        
        
        String outpath = sb.getCurrentUser().getHomeDir() + File.separator + omics1;
        String outpatht = sb.getCurrentUser().getHomeDir() + File.separator + omics2;
        
        
        
        DataUtils.copyFile(new File(inpath), new File(outpath));
        DataUtils.copyFile(new File(inpatht), new File(outpatht));
        
        RDataUtils.readOmicsTable(sb.getRConnection(), omics1, omics2); //do not have to check the error unlike new data upload cuz i know this sample works fine...

        setup16SUploadInfo(omics1, omics2);
    }  
    
     //It used to be doDefaultMammalianAnalysis
    private void doDefaultDILGOMAnalysis() {
        sb.setModuleType("omics");
        sb.setDataFormat("text");
        

        String inpath = ab.getTestrnaDataPath();
        String inpatht = ab.getTestmetaboliteDataPath();
        
        
        String omics1 = DataUtils.getJustFileName(inpath);
        String omics2 = DataUtils.getJustFileName(inpatht);
        
        
        String outpath = sb.getCurrentUser().getHomeDir() + File.separator + omics1;
        String outpatht = sb.getCurrentUser().getHomeDir() + File.separator + omics2;
        
        
        
        DataUtils.copyFile(new File(inpath), new File(outpath));
        DataUtils.copyFile(new File(inpatht), new File(outpatht));
        
        RDataUtils.readOmicsTable(sb.getRConnection(), omics1, omics2);

        setup16SUploadInfo(omics1, omics2);
    }  
    
         //It used to be doDefaultMammalianAnalysis
    private void doDefaultCANCERAnalysis() {
        sb.setModuleType("omics");
        sb.setDataFormat("text");
        

        String inpath = ab.getTestmirnaDataPath();
        String inpatht = ab.getTestmrna2DataPath();
        
        
        String omics1 = DataUtils.getJustFileName(inpath);
        String omics2 = DataUtils.getJustFileName(inpatht);
        
        
        String outpath = sb.getCurrentUser().getHomeDir() + File.separator + omics1;
        String outpatht = sb.getCurrentUser().getHomeDir() + File.separator + omics2;
        
        
        
        DataUtils.copyFile(new File(inpath), new File(outpath));
        DataUtils.copyFile(new File(inpatht), new File(outpatht));
        
        RDataUtils.readOmicsTable(sb.getRConnection(), omics1, omics2);

        setup16SUploadInfo2(omics1, omics2); //no variable sanity check go through before procmsg reading R
    }  
    

    //It used to be doDefaultMammalianAnalysis
    private void doDefaultNUTRIMICEAnalysis() {
        sb.setModuleType("omics");
        sb.setDataFormat("text");
        

        String inpath = ab.getTestlipidDataPath();
        String inpatht = ab.getTestnutrigeneDataPath();
        
        
        String omics1 = DataUtils.getJustFileName(inpath);
        String omics2 = DataUtils.getJustFileName(inpatht);
        
        
        String outpath = sb.getCurrentUser().getHomeDir() + File.separator + omics1;
        String outpatht = sb.getCurrentUser().getHomeDir() + File.separator + omics2;
        
        
        
        DataUtils.copyFile(new File(inpath), new File(outpath));
        DataUtils.copyFile(new File(inpatht), new File(outpatht));
        
        RDataUtils.readOmicsTable(sb.getRConnection(), omics1, omics2);

        setup16SUploadInfo2(omics1, omics2); //no variable sanity check go through before procmsg reading R
    }   
    
    private void doDefaultIBDHUMANAnalysis(){
        sb.setModuleType("omics");
        sb.setDataFormat("text");
        

        String inpath = ab.getTestmetabDataPath();
        String inpatht = ab.getTestmicrobioDataPath();
        
        
        String omics1 = DataUtils.getJustFileName(inpath);
        String omics2 = DataUtils.getJustFileName(inpatht);
        
        
        String outpath = sb.getCurrentUser().getHomeDir() + File.separator + omics1;
        String outpatht = sb.getCurrentUser().getHomeDir() + File.separator + omics2;
        
        
        DataUtils.copyFile(new File(inpath), new File(outpath));
        DataUtils.copyFile(new File(inpatht), new File(outpatht));
        
        RDataUtils.readOmicsTable(sb.getRConnection(), omics1, omics2);

        setup16SUploadInfo2(omics1, omics2); 
    } 
  
}
 


    






