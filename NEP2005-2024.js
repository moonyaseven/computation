/************************************************************
 * 1. 基本配置与省份列表
 ************************************************************/
var startYear = 2005;
var endYear = 2024;

var ProvinceList = ['Aomen','Xianggang','Taiwan',
  'Anhui','Beijing','Chongqing','Fujian','Gansu','Guangdong',
  'Guangxi','Guizhou','Hainan','Hebei','Heilongjiang','Henan',
  'Hubei','Hunan','Jiangsu','Jiangxi','Jilin','Liaoning',
  'Neimenggu','Ningxia','Qinghai','Shan_xi','Shandong','Shanghai',
  'Shanxi','Sichuan','Tianjin', 'Xinjiang','Xizang','Yunnan','Zhejiang'
];

/************************************************************
 * 2. 工具函数定义 
 ************************************************************/
function getDaysInMonth(year, month) {
  var startDate = ee.Date.fromYMD(year, month, 1);
  var endDate = startDate.advance(1, 'month');
  return endDate.difference(startDate, 'day');
}

function calculate_FPAR(ndvi, landcover) {
  var FPAR_min = 0.001;
  var FPAR_max = 0.95;
  var FPAR_delta = FPAR_max - FPAR_min;
  var NDVI_min = ee.Image(0.05);

  var NDVI_max = ee.Image(0.604);
  NDVI_max = NDVI_max.where(landcover.eq(1), 0.604)
                     .where(landcover.eq(2), 0.70)
                     .where(landcover.eq(3), 0.639);

  var NDVI_delta = NDVI_max.subtract(NDVI_min);
  var fpar_ndvi = ndvi.subtract(NDVI_min).divide(NDVI_delta).multiply(FPAR_delta).add(FPAR_min);
  var sr = ee.Image(1).add(ndvi).divide(ee.Image(1).subtract(ndvi));
  var sr_min = ee.Image(1).add(NDVI_min).divide(ee.Image(1).subtract(NDVI_min));
  var sr_max = ee.Image(1).add(NDVI_max).divide(ee.Image(1).subtract(NDVI_max));
  var fpar_sr = sr.subtract(sr_min).divide(sr_max.subtract(sr_min)).multiply(FPAR_delta).add(FPAR_min);

  return fpar_ndvi.add(fpar_sr).multiply(0.5).clamp(0.001, 0.95);
}

function calculate_LUE(t, et, pet, landCover, year, t_opt) {
  var T = t;
  var Tε1 = ee.Image(0.8).add(t_opt.multiply(0.02)).add(t_opt.multiply(t_opt).multiply(-0.0005));
  Tε1 = Tε1.where(T.lte(-10), 0);
  
  var denom1 = t_opt.subtract(10).subtract(T).multiply(0.2).exp().add(1);
  var denom2 = T.subtract(t_opt).subtract(10).multiply(0.3).exp().add(1);
  var Tε2_original = ee.Image(1.184).divide(denom1).divide(denom2);
  
  var extreme = T.gte(t_opt.add(10)).or(T.lte(t_opt.subtract(13)));
  var Tε2 = Tε2_original.where(extreme, Tε2_original.multiply(0.5));
  
  var Wε = et.divide(pet.where(pet.eq(0), 0.01)).multiply(0.5).add(0.5).clamp(0, 1);
  
  var ε_max = ee.Image(0.542);
  ε_max = ε_max.where(landCover.eq(1), 0.542).where(landCover.eq(2), 0.72)
               .where(landCover.eq(3), 0.429).where(landCover.eq(4), 0.542)
               .where(landCover.eq(5), 0).where(landCover.eq(6), 0)
               .where(landCover.eq(7), 0.542).where(landCover.eq(8), 0)
               .where(landCover.eq(9), 0.542);
               
  return Tε1.multiply(Tε2).multiply(Wε).multiply(ε_max);
}

/************************************************************
 * 3. 分省与年份嵌套循环计算
 ************************************************************/
for (var year = startYear; year <= endYear; year++) {
  
  for (var j = 0; j < ProvinceList.length; j++) {
    var loc = ProvinceList[j];
    
    // 获取行政边界 (参照侵蚀脚本逻辑)
    var table;
    if (loc == 'Jilin') {
      table = ee.FeatureCollection("users/moonya20001207/area/Jilin");
    } else {
      table = ee.FeatureCollection("users/moonya20001207/area/" + loc);
    }
    
    print('Current process: ' + loc + ' ' + year);

    var annualNEP = null;
    var scale = 30;

    // 逐月循环计算
    for (var i = 1; i <= 12; i++) {
      var startdate = ee.Date.fromYMD(year, i, 1);
      var enddate = startdate.advance(1, 'month');
      var enddate2 = ee.Date.fromYMD(year, 12, 31);

      // 1. NDVI (Landsat + MODIS 融合)
      var monthlyNDVI_modis = ee.ImageCollection('MODIS/061/MOD13A1')
        .filterDate(startdate, enddate).select('NDVI').max().multiply(0.0001).clip(table)
        .reproject({crs: 'EPSG:4326', scale: scale}).unmask(0);
        
      var ndvi = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_8DAY_NDVI')
        .filterDate(startdate, enddate).select('NDVI').max().clip(table)
        .unmask(monthlyNDVI_modis);


      // 2. ET & PET (支持2023年后自动切换至V6.1)
      function getMonthlyETPET(year, start, end, band, defVal) {
        var collName = (year >= 2023) ? 'MODIS/061/MOD16A2' : 'MODIS/006/MOD16A2';
        var coll = ee.ImageCollection(collName).filterDate(start, end).select(band);
        var hasData = coll.size().gt(0);
        var meanImg = coll.mean().multiply(0.1);
        return ee.Image(ee.Algorithms.If(hasData, meanImg, ee.Image.constant(defVal)));
      }
      function getAnnualETPET(year, band, defVal) {
        var start = ee.Date.fromYMD(year, 1, 1);
        var end = ee.Date.fromYMD(year, 12, 31);
        return getMonthlyETPET(year, start, end, band, defVal);
      }
      
      var et_year = getAnnualETPET(year, 'ET', 0);
      var pet_year = getAnnualETPET(year, 'PET', 0.01);
      var et_month = getMonthlyETPET(year, startdate, enddate, 'ET', 0);
      var pet_month = getMonthlyETPET(year, startdate, enddate, 'PET', 0.01);
      
      var et = et_month.clip(table).unmask(et_year.clip(table));
      var pet = pet_month.clip(table).unmask(pet_year.clip(table));
      
      var et_re = et.resample('bilinear').reproject({crs: "EPSG:4326", scale: scale});
      var pet_re = pet.resample('bilinear').reproject({crs: "EPSG:4326", scale: scale});

      // 3. 气象数据 (SRAD, T, P)
      var days = getDaysInMonth(year, i);
      var srad = ee.Image('users/moonya20001207/SRAD/SRAD_' + year + '_' + i).toFloat().clip(table)
        .multiply(864).divide(1000000).multiply(days).resample('bilinear').reproject({crs: "EPSG:4326", scale: scale});

      var t = ee.Image('users/moonya20001207/temp/temp_' + year + '_' + i).multiply(0.1).clip(table)
        .resample('bilinear').reproject({crs: "EPSG:4326", scale: scale});

      var p = ee.Image('users/moonya20001207/precipitation/Pr_' + year + '_' + i).multiply(0.1).clip(table)
        .resample('bilinear').reproject({crs: "EPSG:4326", scale: scale});

      var t_opt = ee.Image('users/moonya20001207/temp/temp_' + year + '_6')
        .add(ee.Image('users/moonya20001207/temp/temp_' + year + '_7'))
        .add(ee.Image('users/moonya20001207/temp/temp_' + year + '_5'))
        .divide(3).multiply(0.1).clip(table);

      // 4. Landcover & NPP
      var landcover = ee.Image('users/moonya20001207/CLCD/CLCD_v01_' + year + '_albert').clip(table)
        .reproject({crs: "EPSG:4326", scale: scale});
      
      var fpar = calculate_FPAR(ndvi, landcover);
      var eps = calculate_LUE(t, et_re, pet_re, landcover, year, t_opt);
      var npp_month = srad.multiply(fpar).multiply(eps).multiply(0.5);

      // 5. Rh & NEP (逐月减法)
      var LUC = ee.Image(1).where(landcover.eq(5).or(landcover.eq(6)).or(landcover.eq(8)), 0);
      var Rh_month = ee.Image(0.22).multiply((t.multiply(0.0913).exp()).add((p.multiply(0.3145).add(1)).log()))
        .multiply(days).multiply(0.465).multiply(LUC);
      
      var nep_month = npp_month.subtract(Rh_month);

      if (annualNEP == null) { annualNEP = nep_month; } else { annualNEP = annualNEP.add(nep_month); }
    }

    /************************************************************
     * 4. 导出逻辑 
     ************************************************************/
    Export.image.toDrive({
      image: annualNEP.clip(table),
      description: 'NEP_' + loc + '_' + year,
      fileNamePrefix: 'NEP_' + loc + '_' + year,
      folder: 'NEP_' + year + '_30M',
      scale: scale,
      region: table,
      crs: "EPSG:4326",
      maxPixels: 1e13
    });
    
    print(loc + ' ' + year + ' export task started.');
  }
}

print('All provinces NEP calculation finished! Check Tasks panel.');