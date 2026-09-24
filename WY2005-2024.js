/************************************************************
 * WY 产水模型：30米分辨率分省、多年份循环计算与导出
 ************************************************************/

/**************** 1. 年份与省份配置 ******************/
var start = 2005;
var end   = 2024;

var ProvinceList = [
  'Aomen','Xianggang','Taiwan','Anhui','Beijing','Chongqing',
  'Fujian','Gansu','Guangdong','Guangxi','Guizhou','Hainan',
  'Hebei','Heilongjiang','Henan','Hubei','Hunan','Jiangsu',
  'Jiangxi','Jilin','Liaoning','Neimenggu','Ningxia',
  'Qinghai','Shan_xi','Shandong','Shanghai','Shanxi','Sichuan',
  'Tianjin', 'Xinjiang','Xizang','Yunnan','Zhejiang'
];

var scale = 30;
var proj = ee.Projection('EPSG:4326').atScale(scale);


/**************** 2. 基础气候数据集与计算函数 ******************/
var climate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE');

// PAWC 计算函数
function calculate_PAWC(m_s, m_silt, m_c, m_om) {
  return ee.Image(54.509)
    .subtract(m_s.multiply(0.132)).subtract(m_s.pow(2).multiply(0.003))
    .subtract(m_silt.multiply(0.055)).subtract(m_silt.pow(2).multiply(0.006))
    .subtract(m_c.multiply(0.738)).add(m_c.pow(2).multiply(0.007))
    .subtract(m_om.multiply(2.688)).add(m_om.pow(2).multiply(0.501))
    .divide(100).clamp(0, 1);
}

// AWC 计算函数
function calculate_AWC(root_depth, rest_depth, PAWC) {
  return ee.Image(root_depth).min(rest_depth).multiply(PAWC);
}

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation');

// 降雨日判定阈值：日降水 > 0.1 mm
var rainThreshold = 0.1;

// 根据年份计算 N 和 Z
function calculate_Z_02N(year) {

  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = ee.Date.fromYMD(year + 1, 1, 1);

  var rainyDays = chirps
    .filterDate(startDate, endDate)
    .map(function(img) {
      return img.gt(rainThreshold)
        .rename('rain_day')
        .toFloat(); })
    .sum()
    .rename('N');

  var Z = rainyDays
    .multiply(0.2).clamp(1, 30) .rename('Z').toFloat(); 
     return ee.Image.cat([rainyDays, Z]);}

/**************** 3. 循环计算与提交导出 ******************/
for(var year = start; year <= end; year++){
  
  for(var j = 0; j < ProvinceList.length; j++){
    
    var loc = ProvinceList[j];
    
    // 加载省份边界
    var table;
    if(loc == 'Jilin') {
      table = ee.FeatureCollection("users/moonya20001207/area/Jilin");
    } else {
      table = ee.FeatureCollection("users/moonya20001207/area/" + loc);
    }
    
    print('Processing: ' + loc + ' - ' + year);
    
    var startdate = ee.Date.fromYMD(year, 1, 1);
    var enddate   = ee.Date.fromYMD(year + 1, 1, 1);
    
    var rainInfoYear = calculate_Z_02N(year);
    var N_year = rainInfoYear.select('N');
    var Z_year = rainInfoYear.select('Z');

    // ------------------ (1) 土壤数据预处理 (90m重采样至30m) ------------------
    var m_s = ee.Image('projects/ndvi-30meter/assets/sand_0-5cm_90m')
      .divide(100).clip(table).resample('bilinear').reproject(proj);
      
    var m_silt = ee.Image('projects/ndvi-30meter/assets/silt_0-5cm_90m')
      .divide(100).clip(table).resample('bilinear').reproject(proj);
      
    var m_c = ee.Image('projects/ndvi-30meter/assets/clay_0-5cm_90m')
      .divide(100).clip(table).resample('bilinear').reproject(proj);
      
    var m_oc = ee.Image('projects/ndvi-30meter/assets/OC_0-5cm_90m')
      .divide(100).clip(table).resample('bilinear').reproject(proj);

    var OM = m_oc.multiply(1.724);
    var PAWC = calculate_PAWC(m_s, m_silt, m_c, OM).toFloat().rename('PAWC');

    // ------------------ (2) 限制性土层厚度 (250m重采样至30m) ------------------
    var rest_depth = ee.Image('users/moonya20001207/BDRICM_M_250m_ll')
      .clip(table).multiply(10).toFloat().resample('bilinear').reproject(proj);

    // ------------------ (3) CLCD土地利用数据 (30m直接投影) ------------------
    var landCover = ee.Image('users/moonya20001207/CLCD/CLCD_v01_' + year + '_albert')
      .clip(table).toInt().reproject(proj);

    // 参数映射 (利用 CLCD 原始分辨率)
    var kc = landCover.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9], 
      [0.82, 0.82, 0.47, 0.48, 1.00, 0.31, 0.50, 0.30, 1.20]
    ).toFloat().rename('kc');

    var root_depth = landCover.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9],
      [400, 1000, 2000, 300, 1, 1, 1, 1, 200]
    ).toFloat().rename('root_depth');

    var is_veg = landCover.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9],
      [1, 1, 1, 1, 0, 0, 0, 0, 1]
    ).rename('is_veg');

    var AWC = calculate_AWC(root_depth, rest_depth, PAWC).toFloat().rename('AWC');

    // ------------------ (4) 降水数据合成 ------------------
    var P = ee.Image('users/moonya20001207/precipitation/Pr_' + year + '_1');
    for (var m = 2; m <= 12; m++) {
      P = P.add(ee.Image('users/moonya20001207/precipitation/Pr_' + year + '_' + m));
    }
    P = P.clip(table).multiply(0.1).resample('bilinear').reproject(proj).rename('P');
    var P_safe = P.max(1.0);

    // ------------------ (5) 潜在蒸散发 ET0 ------------------
    var ET0 = climate.filterDate(startdate, enddate)
      .select('pet').sum().multiply(0.1)
      .clip(table).toFloat().resample('bilinear').reproject(proj).rename('ET0');

    // ------------------ (6) Budyko 产水逻辑计算 ------------------
    var PET = kc.multiply(ET0).toFloat().rename('PET');
    
     // 年降雨日数 N
    var N = N_year.clip(table).toFloat() .resample('bilinear') .reproject(proj) .rename('N');
    
    // Z = 0.2N
    var Z = Z_year.clip(table).toFloat().resample('bilinear').reproject(proj).rename('Z');
  
    var w = AWC.divide(P_safe).multiply(Z).add(1.25).clamp(1.25, 5).rename('w');
    var pet_p = PET.divide(P_safe).rename('PET_P');

    // 植被与非植被 AET/P 区分计算
    var aet_p_veg = ee.Image(1)
      .add(pet_p)
      .subtract((ee.Image(1).add(pet_p.pow(w))).pow(ee.Image(1).divide(w)))
      .clamp(0, 1);

    var aet_p_nonveg = pet_p.min(1.0);
    var aet_p_ratio = aet_p_nonveg.where(is_veg.eq(1), aet_p_veg).rename('AET_P');

    // 模型估算产水深度 (mm)
    var WY_model = P.multiply(ee.Image(1).subtract(aet_p_ratio))
      .max(0).rename('WY').toFloat();

    // ------------------ (7) 影像导出设置 (与土壤保持保持一致) ------------------
    Export.image.toDrive({
      image: WY_model.clip(table),
      description: 'WY_' + loc + '_' + year,
      fileNamePrefix: 'WY_' + loc + '_' + year,
      folder: 'WY_' + year + '_30M',
      scale: scale,
      region: table,
      crs: "EPSG:4326",
      maxPixels: 1e13
    });

    print(loc + ' - Task generated at 30m resolution.');
  }
}

print('All tasks for Water Yield have been created. Please visit the GEE Tasks tab to run them.');