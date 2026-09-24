/**************** year config ******************/
var start = 2005;
var end = 2024;
            
/**************** 90 M soil data ******************/
var sand = ee.Image('users/dsglsheng/SOIL_90m/btsnd05_90m').divide(10)
var silt = ee.Image('users/dsglsheng/SOIL_90m/btslt05_90m').divide(10)
var clay = ee.Image('users/dsglsheng/SOIL_90m/btcly05_90m').divide(10)
var soc = ee.Image('users/dsglsheng/SOIL_90m/soc05').divide(1000)
            

// SRTM 30m
var dataset = ee.Image('USGS/SRTMGL1_003');
var srtm_elevation = dataset.select('elevation');
var srtm_slope = ee.Terrain.slope(srtm_elevation).rename('slope');


/****************** Rain DownScale source Data Import **********************/
var climate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')

var aet = climate.select('aet')
var def = climate.select('def')
var pdsi = climate.select('pdsi')
var ro = climate.select('ro')
var soil = climate.select('soil')
var srad = climate.select('srad')
var swe = climate.select('swe')
var tmmn = climate.select('tmmn')
var tmmx = climate.select('tmmx')
var vap = climate.select('vap')
var vpd = climate.select('vpd')
var vs = climate.select('vs')

var dem = ee.Image('CGIAR/SRTM90_V4').select('elevation');

var slope = ee.Terrain.slope(dem).rename('slope');
/****************** Rain DownScale source Data Import Complete! **********************/


/**************** mosaic the LS dataset ******************/
var list = ee.List([])
for(var i = 1; i < 9; i++){list = list.add(ee.Image('users/moonya20001207/LS/ls_' + i))};
var col = ee.ImageCollection.fromImages(list)
var merge = col.mosaic()


/**************** split the province ******************/
var ProvinceList = ['Aomen','Xianggang','Taiwan','Anhui','Beijing','Chongqing',
'Fujian','Gansu','Guangdong','Guangxi','Guizhou','Hainan',
'Hebei','Heilongjiang','Henan','Hubei','Hunan','Jiangsu',
'Jiangxi','Jilin','Liaoning','Neimenggu','Ningxia',
'Qinghai','Shan_xi','Shandong','Shanghai','Shanxi','Sichuan',
'Tianjin', 'Xinjiang','Xizang','Yunnan',
'Zhejiang']
// var ProvinceList = ['China_shp']



/*****************************************************************
 * process Landsat images
*****************************************************************/
// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              // .addBands(thermalBands, null, true);
}

// reomove cloud for Landsat-8
function cloudRemoval(image) { 
  var cloudShadowBitMask = (1 << 4); 
  var cloudsBitMask = (1 << 3); 
  var qa = image.select('QA_PIXEL'); 
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0) 
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0)); 
  var mask2 = image.select("blue").gt(0.2);                 
  return image.updateMask(mask).updateMask(mask2.not()).toDouble()
              .copyProperties(image)
              .copyProperties(image, ["system:time_start"]);
} 


// Assign a common name to the sensor-specific bands.
var LC9_BANDS = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL']; //Landsat 9
var LC8_BANDS = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL']; //Landsat 8
var LC7_BANDS = ['SR_B1','SR_B2','SR_B3', 'SR_B4','SR_B5','SR_B7','QA_PIXEL']; //Landsat 7
var LC5_BANDS = ['SR_B1','SR_B2','SR_B3', 'SR_B4','SR_B5','SR_B7','QA_PIXEL']; //Landsat 5
var S2_BANDS  = ['B2','B3','B4','B8','B11','B12']; // Sentinel-2
var STD_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2','QA_PIXEL'];


/*****************************************************************
 * calculate the R factor
*****************************************************************/
function cal_R(rain_image) {
  var r_factor = (rain_image.pow(1.6548)).multiply(0.0534)
  return r_factor.rename('R_Factor')
}


/*****************************************************************
 * calculate the K factor
*****************************************************************/
function cal_K(SAN, SIL, CLA, C, SN1) {
  var A = ((((SAN.multiply((SIL.divide(-100)).add(1))).multiply(-0.0256)).exp()).multiply(0.3)).add(0.2)
  var B = SIL.divide(CLA.add(SIL))
  var C_ = (C.multiply(0.25)).divide(C.add(((C.multiply(-2.95)).add(3.72)).exp()))
  var D = (SN1.multiply(0.7)).divide(SN1.add(((SN1.multiply(22.9)).add(-5.51)).exp()))
  var k_factor = SAN.expression("A*(B**0.3)*(1-C)*(1-D)", {
      "A": A,
      "B": B,
      "C": C_,
      "D": D,
  })
  // .rename('K_Factor').multiply(0.1317);
  k_factor = ((k_factor.rename('K_Factor').multiply(0.51575)).add(-0.01383)).multiply(0.1317)
  return k_factor;
}

/*****************************************************************
 * calculate the FVC and NDVI
*****************************************************************/
//*****NDVI & FVC*****
function FVC(img){
  var ndvi = img.select("NDVI");
 var fvc = img.expression(
   "108.49 * NDVI + 0.717",
   {
     "NDVI": ndvi,
   }
 ).rename('FVC').divide(100);
 return fvc;
}

function NDVI(img) {
 var nir = img.select("nir");
 var red = img.select("red");
 var ndvi = img.expression(
   "(nir - red)/(nir + red)",
   {
     "nir": nir,
     "red": red
   }
 ).rename("NDVI");
 return ndvi;
}


function cal_C_Factor(landCover, FVC, NDVI) {
  var C = landCover.expression(
    'land_cover == 1 ? Cropland_cal : ' +
    'land_cover == 2 ? Forest_cal : ' +
    'land_cover == 3 ? Shrub_cal : ' +
    'land_cover == 4 ? Grass_cal : ' +
    'land_cover == 7 ? Barren_cal : ' +
    'default_formula',
    {
      'land_cover': landCover,
      // 'Cropland_cal': FVC.expression('1 - fvc', {'fvc': FVC.select('FVC')}),
      'Cropland_cal': FVC.expression('1 - fvc', {'fvc': FVC.select('FVC')}),
      'Forest_cal': FVC.expression(
        '(fvc < 0.1) ? 0.1 : ' +
        '(fvc < 0.3) ? 0.08 : ' +
        '(fvc < 0.5) ? 0.06 : ' +
        '(fvc < 0.7) ? 0.02 : ' +
        '(fvc < 0.9) ? 0.004 : ' +
        '0.001',
        {
          'fvc': FVC.select('FVC'),
          'ndvi': NDVI.select('NDVI')
        }
      ),
      'Shrub_cal': FVC.expression('min+(max - min) * (1 - fvc)', {'min': 0.01, 'max': 0.15, 'fvc': FVC.select('FVC')}),
      'Grass_cal': FVC.expression('(1 - ndvi) / 2', {'ndvi': NDVI.select('NDVI')}),
      'Barren_cal': FVC.expression('min+(max - min) * (1 - fvc)', {'min': 0.1, 'max': 0.5, 'fvc': FVC.select('FVC')}),
      'default_formula': ee.Image(0)
    }
  );
  
  return C.rename('C_Factor');
}


/*****************************************************************
        * provincal loop calculate
*****************************************************************/

for(var year=start; year<=end; year++){
  var RList = ee.List([])
  
  for(var j=0;j<ProvinceList.length;j++){
    
    var loc = ProvinceList[j]
    
    if(loc == 'Jilin') {
      print("Jilin come on!")
      var table = ee.FeatureCollection("users/moonya20001207/area/Jilin");
    } else {
      var table = ee.FeatureCollection("users/moonya20001207/area/"+loc);
    }
    
    print('Current process is '+ loc + year)
    var landCover = ee.Image('users/moonya20001207/CLCD/CLCD_v01_'+ year +'_albert').clip(table)
      
    
    /*****************************************************************
        * C factor
    *****************************************************************/
    var L9Col = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                    .filterBounds(table)
                    .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                    .map(applyScaleFactors)
                    .select(LC9_BANDS,STD_NAMES)
                    .map(cloudRemoval);
    var L8Col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                    .filterBounds(table)
                    .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                    .map(applyScaleFactors)
                    .select(LC8_BANDS,STD_NAMES)
                    .map(cloudRemoval);
    var L7Col = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                    .filterBounds(table)
                    .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                    .map(applyScaleFactors)
                    .select(LC7_BANDS,STD_NAMES)
                    .map(cloudRemoval);
    var L5Col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                    .filterBounds(table)
                    .filter(ee.Filter.date(year+'-01-01', year+'-12-31'))
                    .map(applyScaleFactors)
                    .select(LC5_BANDS,STD_NAMES)
                    .map(cloudRemoval);
                    
    var LandsatCol = ee.ImageCollection(L9Col.merge(L8Col).merge(L7Col).merge(L5Col))
    
    /** caculate the NDVI **/
    var final_image = LandsatCol.mean().clip(table);
    var image_ndvi = NDVI(final_image)
    var ndvi = image_ndvi.select('NDVI').clip(table)
    var ndvi = ndvi.where(ndvi.lt(-1), -1);
    var ndvi = ndvi.where(ndvi.gt(1), 1);
    var fvc = FVC(ndvi).clip(table)
    
    var C_Factor = cal_C_Factor(landCover, fvc, ndvi).clip(table)
    
    
    /*****************************************************************
        * R factor
    *****************************************************************/
    var scale = 30;
    var P = ee.Image('users/moonya20001207/precipitation/Pr_' + year + '_1');
    for (var m = 2; m <= 12; m++) {
      P = P.add(ee.Image('users/moonya20001207/precipitation/Pr_' + year + '_' + m));
    }
    P = P.clip(table).multiply(0.1).reproject({crs:"EPSG:4326", scale:30 }).rename('P');
    var rain_image = P.max(1.0); 
    var R_Factor = cal_R(rain_image).clip(table);
    
    
    /*****************************************************************
        * K factor
    *****************************************************************/
    var SAN = sand.clip(table)
    var SIL = silt.clip(table)
    var CLA = clay.clip(table)
    var C = soc.clip(table)
    var SN1 = (SAN.multiply(-1).add(1)).divide(100)
    
    var K_Factor = cal_K(SAN, SIL, CLA, C, SN1).clip(table)
    
    if(scale < 90){
      var K_Factor = K_Factor.resample('bilinear').reproject({
        crs: K_Factor.projection(),
        scale: scale
      });
    }
    K_Factor = K_Factor.clip(table)
    
    /*****************************************************************
        * LS factor
    *****************************************************************/
    var LS_Factor = merge.clip(table)
    
    /*****************************************************************
        * P factor
    *****************************************************************/
    var local_slope = srtm_slope.clip(table)
    
    landCover = landCover.rename('lulc')
    local_slope = local_slope.rename('slope')
    var lulc_slope = landCover.addBands(local_slope)
    
    var P_Factor = lulc_slope.expression(

      "(b('lulc') == 2) ? 1" +
      ": (b('lulc') == 7) ? 1" +
      ": (b('lulc') == 3) ? 0.29" +
      ": (b('lulc') == 4) ? 0.41" +
      ": (b('slope') < 5) and (b('lulc')==1) ? 0.49" +
      ": (b('slope') < 7) and (b('lulc')==1) ? 0.59" +
      ": (b('slope') < 9) and (b('lulc')==1) ? 0.65" +
      ": (b('slope') < 12) and (b('lulc')==1) ? 0.70" +
      ": (b('slope') < 20) and (b('lulc')==1) ? 0.81" +
      ": (b('slope') < 24) and (b('lulc')==1) ? 0.95" +
      ": (b('slope') > 24) and (b('lulc')==1) ? 1" +
      ": 0"
    ).rename('P').clip(table);
    
  
    /*****************************************************************
        * Soil Erosion A
    *****************************************************************/
  var ASE = C_Factor.multiply(R_Factor).multiply(LS_Factor).multiply(K_Factor).multiply(P_Factor)
    .clip(table).rename('ASE');
  var PSE = R_Factor.multiply(LS_Factor).multiply(K_Factor).clip(table).rename('PSE');
  var SR = PSE.subtract(ASE).rename('SR');
  var SR_1 = SR.where(SR.lt(0), 0).where(SR.gte(1e9), 0).clip(table);
    
    
    // var scale = 30
    
    // Export.image.toDrive({
    //   image: C_Factor.clip(table),
    //   description: 'C_'+loc+'_'+year,
    //   fileNamePrefix: 'C_'+loc+'_'+year,
    //   folder:'C_Factor',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    // Export.image.toDrive({
    //   image: K_Factor.clip(table),
    //   description: 'K_'+loc,
    //   fileNamePrefix: 'K_'+loc,
    //   folder:'K_Factor',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    // Export.image.toDrive({
    //   image: R_Factor.clip(table),
    //   description: 'R_'+loc+'_'+year,
    //   fileNamePrefix: 'R_'+loc+'_'+year,
    //   folder:'R_Factor',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    // Export.image.toDrive({
    //   image: LS_Factor.clip(table),
    //   description: 'LS_'+loc,
    //   fileNamePrefix: 'LS_'+loc,
    //   folder:'LS_Factor_1KM',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    
    // Export.image.toDrive({
    //   image: P_Factor.clip(table),
    //   description: 'P_'+loc+'_'+year,
    //   fileNamePrefix: 'P_'+loc+'_'+year,
    //   folder:'P_Factor_1KM',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    // Export.image.toDrive({
    //   image: A_1.clip(table),
    //   description: 'SE_'+loc+'_'+year,
    //   fileNamePrefix: 'SE_'+loc+'_'+year,
    //   folder:'SE_' + year + '_1KM',
    //   scale: scale,
    //   region: table,
    //   crs: "EPSG:4326",
    //   maxPixels:1e13
    // });
    
    // var pixel = A_1.sampleRegions({
    //   collection: point,
    //   scale: scale,
    //   geometries: true
    // })
    
    //Export.image.toDrive({
    //  image: SE_1.clip(table),
    // description: 'SE_'+loc+'_'+year,
    //  fileNamePrefix: 'SE_'+loc+'_'+year,
    //  folder:'SE_' + year + '_30M',
    //  scale: scale,
   //   region: table,
   //   crs: "EPSG:4326",
    //  maxPixels:1e13
   // });
    Export.image.toDrive({
      image: SR_1.clip(table),
      description: 'SR_'+loc+'_'+year,
      fileNamePrefix: 'SR_'+loc+'_'+year,
      folder:'SR_' + year + '_30M',
      scale: scale,
      region: table,
      crs: "EPSG:4326",
      maxPixels:1e13
    });
print(loc, scale)
  }
}

print('calculate finish! start export the csv')