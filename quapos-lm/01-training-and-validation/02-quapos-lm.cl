/*
OpenCL RandomForestClassifier
classifier_class_name = ObjectSegmenter
feature_specification = gaussian_blur=1 difference_of_gaussian=1 laplace_box_of_gaussian_blur=1
num_ground_truth_dimensions = 3
num_classes = 2
num_features = 3
max_depth = 2
num_trees = 100
feature_importances = 0.32557488170342097,0.4231073391932076,0.25131777910337144
positive_class_identifier = 2
apoc_version = 0.12.0
*/
__kernel void predict (IMAGE_in0_TYPE in0, IMAGE_in1_TYPE in1, IMAGE_in2_TYPE in2, IMAGE_out_TYPE out) {
 sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
 const int x = get_global_id(0);
 const int y = get_global_id(1);
 const int z = get_global_id(2);
 float i0 = READ_IMAGE(in0, sampler, POS_in0_INSTANCE(x,y,z,0)).x;
 float i1 = READ_IMAGE(in1, sampler, POS_in1_INSTANCE(x,y,z,0)).x;
 float i2 = READ_IMAGE(in2, sampler, POS_in2_INSTANCE(x,y,z,0)).x;
 float s0=0;
 float s1=0;
if(i2<-10111.8388671875){
 if(i0<454.7654724121094){
  s0+=168.0;
  s1+=49.0;
 } else {
  s0+=97.0;
  s1+=1846.0;
 }
} else {
 if(i0<354.4898681640625){
  s0+=1947.0;
  s1+=51.0;
 } else {
  s0+=8.0;
  s1+=117.0;
 }
}
if(i0<406.37890625){
 if(i2<-6836.1669921875){
  s0+=613.0;
  s1+=89.0;
 } else {
  s0+=1406.0;
  s1+=10.0;
 }
} else {
 if(i1<-18.950393676757812){
  s0+=126.0;
  s1+=109.0;
 } else {
  s0+=19.0;
  s1+=1911.0;
 }
}
if(i0<397.7173156738281){
 if(i2<-6836.1669921875){
  s0+=616.0;
  s1+=79.0;
 } else {
  s0+=1402.0;
  s1+=8.0;
 }
} else {
 if(i1<-18.281112670898438){
  s0+=127.0;
  s1+=107.0;
 } else {
  s0+=29.0;
  s1+=1915.0;
 }
}
if(i1<5.102691650390625){
 if(i2<-12339.9755859375){
  s0+=107.0;
  s1+=225.0;
 } else {
  s0+=2049.0;
  s1+=52.0;
 }
} else {
 if(i2<-6058.4921875){
  s0+=17.0;
  s1+=1821.0;
 } else {
  s0+=12.0;
 }
}
if(i0<362.9952392578125){
 if(i2<-6877.66748046875){
  s0+=549.0;
  s1+=54.0;
 } else {
  s0+=1444.0;
  s1+=9.0;
 }
} else {
 if(i0<533.5504150390625){
  s0+=146.0;
  s1+=265.0;
 } else {
  s0+=42.0;
  s1+=1774.0;
 }
}
if(i1<3.3112564086914062){
 if(i2<-13857.537109375){
  s0+=58.0;
  s1+=204.0;
 } else {
  s0+=2114.0;
  s1+=47.0;
 }
} else {
 if(i1<6.5342559814453125){
  s0+=21.0;
  s1+=36.0;
 } else {
  s0+=17.0;
  s1+=1786.0;
 }
}
if(i2<-10154.09765625){
 if(i1<-23.927459716796875){
  s0+=232.0;
  s1+=105.0;
 } else {
  s0+=35.0;
  s1+=1845.0;
 }
} else {
 if(i1<6.9555206298828125){
  s0+=1894.0;
  s1+=16.0;
 } else {
  s0+=8.0;
  s1+=148.0;
 }
}
if(i1<4.6247711181640625){
 if(i1<-2.818115234375){
  s0+=2041.0;
  s1+=182.0;
 } else {
  s0+=127.0;
  s1+=59.0;
 }
} else {
 if(i1<12.165298461914062){
  s0+=37.0;
  s1+=105.0;
 } else {
  s0+=6.0;
  s1+=1726.0;
 }
}
if(i1<5.102691650390625){
 if(i2<-13803.9208984375){
  s0+=78.0;
  s1+=210.0;
 } else {
  s0+=2127.0;
  s1+=56.0;
 }
} else {
 if(i2<-5986.8583984375){
  s0+=17.0;
  s1+=1786.0;
 } else {
  s0+=8.0;
  s1+=1.0;
 }
}
if(i0<356.5364990234375){
 if(i0<322.1090393066406){
  s0+=1871.0;
  s1+=16.0;
 } else {
  s0+=103.0;
  s1+=31.0;
 }
} else {
 if(i0<486.2808837890625){
  s0+=141.0;
  s1+=199.0;
 } else {
  s0+=68.0;
  s1+=1854.0;
 }
}
if(i1<6.43902587890625){
 if(i0<474.216064453125){
  s0+=2078.0;
  s1+=32.0;
 } else {
  s0+=72.0;
  s1+=262.0;
 }
} else {
 if(i0<303.9437561035156){
  s0+=11.0;
  s1+=6.0;
 } else {
  s0+=3.0;
  s1+=1819.0;
 }
}
if(i1<5.154563903808594){
 if(i2<-13803.9208984375){
  s0+=78.0;
  s1+=217.0;
 } else {
  s0+=2089.0;
  s1+=61.0;
 }
} else {
 if(i1<15.837936401367188){
  s0+=32.0;
  s1+=134.0;
 } else {
  s0+=5.0;
  s1+=1667.0;
 }
}
if(i2<-9301.8828125){
 if(i0<429.4349365234375){
  s0+=202.0;
  s1+=47.0;
 } else {
  s0+=114.0;
  s1+=1962.0;
 }
} else {
 if(i0<321.6493225097656){
  s0+=1808.0;
  s1+=8.0;
 } else {
  s0+=42.0;
  s1+=100.0;
 }
}
if(i2<-10116.294921875){
 if(i1<-18.464279174804688){
  s0+=228.0;
  s1+=120.0;
 } else {
  s0+=23.0;
  s1+=1849.0;
 }
} else {
 if(i2<-7278.2216796875){
  s0+=400.0;
  s1+=127.0;
 } else {
  s0+=1516.0;
  s1+=20.0;
 }
}
if(i0<365.0766906738281){
 if(i2<-7350.8330078125){
  s0+=426.0;
  s1+=39.0;
 } else {
  s0+=1595.0;
  s1+=8.0;
 }
} else {
 if(i0<480.5464172363281){
  s0+=122.0;
  s1+=189.0;
 } else {
  s0+=68.0;
  s1+=1836.0;
 }
}
if(i2<-10111.8388671875){
 if(i0<462.20245361328125){
  s0+=169.0;
  s1+=58.0;
 } else {
  s0+=71.0;
  s1+=1938.0;
 }
} else {
 if(i0<353.0996398925781){
  s0+=1890.0;
  s1+=33.0;
 } else {
  s0+=17.0;
  s1+=107.0;
 }
}
if(i1<4.6247711181640625){
 if(i2<-13395.0556640625){
  s0+=78.0;
  s1+=220.0;
 } else {
  s0+=2149.0;
  s1+=62.0;
 }
} else {
 if(i2<-7115.203125){
  s0+=10.0;
  s1+=1737.0;
 } else {
  s0+=19.0;
  s1+=8.0;
 }
}
if(i0<365.1287841796875){
 if(i2<-6174.80126953125){
  s0+=774.0;
  s1+=55.0;
 } else {
  s0+=1296.0;
 }
} else {
 if(i1<-18.41778564453125){
  s0+=142.0;
  s1+=115.0;
 } else {
  s0+=24.0;
  s1+=1877.0;
 }
}
if(i1<3.3112564086914062){
 if(i0<474.216064453125){
  s0+=2068.0;
  s1+=37.0;
 } else {
  s0+=73.0;
  s1+=215.0;
 }
} else {
 if(i0<310.71514892578125){
  s0+=27.0;
  s1+=11.0;
 } else {
  s0+=18.0;
  s1+=1834.0;
 }
}
if(i1<2.80157470703125){
 if(i0<457.3300476074219){
  s0+=2128.0;
  s1+=31.0;
 } else {
  s0+=76.0;
  s1+=200.0;
 }
} else {
 if(i2<-6380.5908203125){
  s0+=19.0;
  s1+=1805.0;
 } else {
  s0+=22.0;
  s1+=2.0;
 }
}
if(i2<-10308.70703125){
 if(i2<-13275.8076171875){
  s0+=83.0;
  s1+=1652.0;
 } else {
  s0+=136.0;
  s1+=259.0;
 }
} else {
 if(i2<-8562.2666015625){
  s0+=212.0;
  s1+=121.0;
 } else {
  s0+=1766.0;
  s1+=54.0;
 }
}
if(i0<364.4610595703125){
 if(i1<7.429473876953125){
  s0+=1946.0;
  s1+=13.0;
 } else {
  s0+=10.0;
  s1+=53.0;
 }
} else {
 if(i1<-18.281112670898438){
  s0+=169.0;
  s1+=98.0;
 } else {
  s0+=27.0;
  s1+=1967.0;
 }
}
if(i0<402.4314880371094){
 if(i0<322.1090393066406){
  s0+=1887.0;
  s1+=8.0;
 } else {
  s0+=177.0;
  s1+=89.0;
 }
} else {
 if(i1<-18.234619140625){
  s0+=129.0;
  s1+=104.0;
 } else {
  s0+=22.0;
  s1+=1867.0;
 }
}
if(i0<397.90911865234375){
 if(i2<-6877.66748046875){
  s0+=588.0;
  s1+=77.0;
 } else {
  s0+=1473.0;
  s1+=10.0;
 }
} else {
 if(i0<533.1624755859375){
  s0+=120.0;
  s1+=282.0;
 } else {
  s0+=49.0;
  s1+=1684.0;
 }
}
if(i1<4.377410888671875){
 if(i2<-14400.20703125){
  s0+=41.0;
  s1+=190.0;
 } else {
  s0+=2141.0;
  s1+=71.0;
 }
} else {
 if(i0<308.5609130859375){
  s0+=13.0;
  s1+=5.0;
 } else {
  s0+=19.0;
  s1+=1803.0;
 }
}
if(i0<404.2063903808594){
 if(i0<322.7826232910156){
  s0+=1861.0;
  s1+=13.0;
 } else {
  s0+=195.0;
  s1+=75.0;
 }
} else {
 if(i0<493.8672790527344){
  s0+=87.0;
  s1+=159.0;
 } else {
  s0+=51.0;
  s1+=1842.0;
 }
}
if(i0<406.37890625){
 if(i1<7.429473876953125){
  s0+=2065.0;
  s1+=12.0;
 } else {
  s0+=15.0;
  s1+=82.0;
 }
} else {
 if(i1<-18.234619140625){
  s0+=117.0;
  s1+=116.0;
 } else {
  s0+=18.0;
  s1+=1858.0;
 }
}
if(i2<-11409.63671875){
 if(i2<-16109.0419921875){
  s0+=32.0;
  s1+=1371.0;
 } else {
  s0+=117.0;
  s1+=379.0;
 }
} else {
 if(i0<362.1824951171875){
  s0+=2039.0;
  s1+=60.0;
 } else {
  s0+=62.0;
  s1+=223.0;
 }
}
if(i2<-10121.10546875){
 if(i2<-12908.51953125){
  s0+=84.0;
  s1+=1633.0;
 } else {
  s0+=165.0;
  s1+=264.0;
 }
} else {
 if(i2<-7278.2216796875){
  s0+=427.0;
  s1+=147.0;
 } else {
  s0+=1555.0;
  s1+=8.0;
 }
}
if(i0<399.22314453125){
 if(i1<7.672821044921875){
  s0+=2076.0;
  s1+=18.0;
 } else {
  s0+=11.0;
  s1+=73.0;
 }
} else {
 if(i0<549.8062744140625){
  s0+=104.0;
  s1+=263.0;
 } else {
  s0+=36.0;
  s1+=1702.0;
 }
}
if(i1<2.8080406188964844){
 if(i2<-12104.111328125){
  s0+=109.0;
  s1+=222.0;
 } else {
  s0+=2016.0;
  s1+=30.0;
 }
} else {
 if(i0<310.8768310546875){
  s0+=30.0;
  s1+=12.0;
 } else {
  s0+=14.0;
  s1+=1850.0;
 }
}
if(i2<-10153.9267578125){
 if(i1<-23.878173828125){
  s0+=196.0;
  s1+=85.0;
 } else {
  s0+=37.0;
  s1+=1894.0;
 }
} else {
 if(i1<7.52520751953125){
  s0+=1897.0;
  s1+=6.0;
 } else {
  s0+=10.0;
  s1+=158.0;
 }
}
if(i1<5.600189208984375){
 if(i0<474.06561279296875){
  s0+=2097.0;
  s1+=42.0;
 } else {
  s0+=82.0;
  s1+=274.0;
 }
} else {
 if(i0<280.22613525390625){
  s0+=10.0;
  s1+=1.0;
 } else {
  s0+=10.0;
  s1+=1767.0;
 }
}
if(i0<402.25018310546875){
 if(i2<-6510.4794921875){
  s0+=730.0;
  s1+=86.0;
 } else {
  s0+=1335.0;
  s1+=3.0;
 }
} else {
 if(i1<-23.927459716796875){
  s0+=99.0;
  s1+=109.0;
 } else {
  s0+=18.0;
  s1+=1903.0;
 }
}
if(i1<5.591796875){
 if(i0<496.076416015625){
  s0+=2069.0;
  s1+=48.0;
 } else {
  s0+=58.0;
  s1+=235.0;
 }
} else {
 if(i0<303.9437561035156){
  s0+=13.0;
  s1+=7.0;
 } else {
  s0+=13.0;
  s1+=1840.0;
 }
}
if(i2<-9823.845703125){
 if(i0<454.9073486328125){
  s0+=189.0;
  s1+=51.0;
 } else {
  s0+=86.0;
  s1+=1901.0;
 }
} else {
 if(i1<6.007179260253906){
  s0+=1895.0;
  s1+=9.0;
 } else {
  s0+=16.0;
  s1+=136.0;
 }
}
if(i0<397.5604248046875){
 if(i2<-6163.68115234375){
  s0+=792.0;
  s1+=87.0;
 } else {
  s0+=1248.0;
  s1+=1.0;
 }
} else {
 if(i1<-23.890045166015625){
  s0+=120.0;
  s1+=84.0;
 } else {
  s0+=33.0;
  s1+=1918.0;
 }
}
if(i1<6.4498443603515625){
 if(i0<493.47802734375){
  s0+=2164.0;
  s1+=48.0;
 } else {
  s0+=57.0;
  s1+=221.0;
 }
} else {
 if(i2<-6020.0703125){
  s0+=7.0;
  s1+=1781.0;
 } else {
  s0+=4.0;
  s1+=1.0;
 }
}
if(i1<2.795604705810547){
 if(i1<-2.8140907287597656){
  s0+=2032.0;
  s1+=191.0;
 } else {
  s0+=119.0;
  s1+=62.0;
 }
} else {
 if(i1<7.688507080078125){
  s0+=27.0;
  s1+=53.0;
 } else {
  s0+=11.0;
  s1+=1788.0;
 }
}
if(i2<-10211.46875){
 if(i2<-12421.400390625){
  s0+=100.0;
  s1+=1731.0;
 } else {
  s0+=133.0;
  s1+=178.0;
 }
} else {
 if(i1<7.429473876953125){
  s0+=1969.0;
  s1+=12.0;
 } else {
  s0+=10.0;
  s1+=150.0;
 }
}
if(i1<2.0205535888671875){
 if(i2<-13882.1396484375){
  s0+=73.0;
  s1+=208.0;
 } else {
  s0+=2099.0;
  s1+=63.0;
 }
} else {
 if(i2<-7602.244140625){
  s0+=20.0;
  s1+=1771.0;
 } else {
  s0+=29.0;
  s1+=20.0;
 }
}
if(i2<-10118.388671875){
 if(i1<-16.044631958007812){
  s0+=208.0;
  s1+=109.0;
 } else {
  s0+=14.0;
  s1+=1777.0;
 }
} else {
 if(i0<354.4898681640625){
  s0+=2006.0;
  s1+=38.0;
 } else {
  s0+=14.0;
  s1+=117.0;
 }
}
if(i1<3.3112564086914062){
 if(i1<-2.280517578125){
  s0+=2031.0;
  s1+=203.0;
 } else {
  s0+=101.0;
  s1+=56.0;
 }
} else {
 if(i1<8.402206420898438){
  s0+=23.0;
  s1+=49.0;
 } else {
  s0+=17.0;
  s1+=1803.0;
 }
}
if(i0<405.945556640625){
 if(i1<7.429473876953125){
  s0+=2102.0;
  s1+=17.0;
 } else {
  s0+=11.0;
  s1+=103.0;
 }
} else {
 if(i2<-16973.927734375){
  s0+=24.0;
  s1+=1257.0;
 } else {
  s0+=120.0;
  s1+=649.0;
 }
}
if(i2<-10156.130859375){
 if(i1<-19.33721923828125){
  s0+=223.0;
  s1+=112.0;
 } else {
  s0+=18.0;
  s1+=1844.0;
 }
} else {
 if(i0<349.28314208984375){
  s0+=1896.0;
  s1+=37.0;
 } else {
  s0+=28.0;
  s1+=125.0;
 }
}
if(i1<6.007179260253906){
 if(i2<-12284.755859375){
  s0+=97.0;
  s1+=235.0;
 } else {
  s0+=2012.0;
  s1+=40.0;
 }
} else {
 if(i1<17.206680297851562){
  s0+=21.0;
  s1+=139.0;
 } else {
  s0+=1.0;
  s1+=1738.0;
 }
}
if(i0<360.99322509765625){
 if(i0<322.65380859375){
  s0+=1881.0;
  s1+=12.0;
 } else {
  s0+=117.0;
  s1+=41.0;
 }
} else {
 if(i2<-18610.484375){
  s0+=11.0;
  s1+=1208.0;
 } else {
  s0+=191.0;
  s1+=822.0;
 }
}
if(i1<4.318756103515625){
 if(i2<-13803.9208984375){
  s0+=55.0;
  s1+=182.0;
 } else {
  s0+=2183.0;
  s1+=57.0;
 }
} else {
 if(i1<6.4498443603515625){
  s0+=15.0;
  s1+=28.0;
 } else {
  s0+=15.0;
  s1+=1748.0;
 }
}
if(i0<399.22314453125){
 if(i2<-6834.109375){
  s0+=579.0;
  s1+=81.0;
 } else {
  s0+=1491.0;
  s1+=5.0;
 }
} else {
 if(i1<-18.281112670898438){
  s0+=145.0;
  s1+=113.0;
 } else {
  s0+=33.0;
  s1+=1836.0;
 }
}
if(i1<6.007179260253906){
 if(i1<-1.9239578247070312){
  s0+=2044.0;
  s1+=200.0;
 } else {
  s0+=131.0;
  s1+=80.0;
 }
} else {
 if(i2<-6131.056640625){
  s0+=22.0;
  s1+=1796.0;
 } else {
  s0+=9.0;
  s1+=1.0;
 }
}
if(i2<-10156.130859375){
 if(i2<-13367.9609375){
  s0+=91.0;
  s1+=1631.0;
 } else {
  s0+=177.0;
  s1+=291.0;
 }
} else {
 if(i0<360.9933776855469){
  s0+=1933.0;
  s1+=60.0;
 } else {
  s0+=4.0;
  s1+=96.0;
 }
}
if(i1<6.061004638671875){
 if(i1<-1.938079833984375){
  s0+=2026.0;
  s1+=217.0;
 } else {
  s0+=121.0;
  s1+=71.0;
 }
} else {
 if(i2<-6131.056640625){
  s0+=17.0;
  s1+=1816.0;
 } else {
  s0+=15.0;
 }
}
if(i1<6.4498443603515625){
 if(i1<-2.376190185546875){
  s0+=2037.0;
  s1+=190.0;
 } else {
  s0+=134.0;
  s1+=85.0;
 }
} else {
 if(i1<17.206680297851562){
  s0+=14.0;
  s1+=139.0;
 } else {
  s1+=1684.0;
 }
}
if(i1<5.102691650390625){
 if(i1<-2.110076904296875){
  s0+=2061.0;
  s1+=199.0;
 } else {
  s0+=122.0;
  s1+=64.0;
 }
} else {
 if(i1<18.222412109375){
  s0+=26.0;
  s1+=169.0;
 } else {
  s0+=1.0;
  s1+=1641.0;
 }
}
if(i1<5.102691650390625){
 if(i1<-2.280517578125){
  s0+=2045.0;
  s1+=231.0;
 } else {
  s0+=120.0;
  s1+=66.0;
 }
} else {
 if(i0<315.7391357421875){
  s0+=20.0;
  s1+=10.0;
 } else {
  s0+=14.0;
  s1+=1777.0;
 }
}
if(i0<397.7173156738281){
 if(i1<7.672821044921875){
  s0+=2053.0;
  s1+=23.0;
 } else {
  s0+=12.0;
  s1+=55.0;
 }
} else {
 if(i1<-18.281112670898438){
  s0+=131.0;
  s1+=124.0;
 } else {
  s0+=13.0;
  s1+=1872.0;
 }
}
if(i1<3.3112564086914062){
 if(i1<-4.3166656494140625){
  s0+=2010.0;
  s1+=169.0;
 } else {
  s0+=150.0;
  s1+=72.0;
 }
} else {
 if(i1<7.672821044921875){
  s0+=35.0;
  s1+=56.0;
 } else {
  s0+=15.0;
  s1+=1776.0;
 }
}
if(i2<-10156.7958984375){
 if(i2<-14839.46875){
  s0+=49.0;
  s1+=1465.0;
 } else {
  s0+=175.0;
  s1+=445.0;
 }
} else {
 if(i1<6.514984130859375){
  s0+=1957.0;
  s1+=18.0;
 } else {
  s0+=15.0;
  s1+=159.0;
 }
}
if(i2<-9764.26953125){
 if(i2<-13367.9609375){
  s0+=72.0;
  s1+=1676.0;
 } else {
  s0+=171.0;
  s1+=291.0;
 }
} else {
 if(i2<-8005.505859375){
  s0+=238.0;
  s1+=100.0;
 } else {
  s0+=1701.0;
  s1+=34.0;
 }
}
if(i1<2.7891387939453125){
 if(i2<-12206.0703125){
  s0+=106.0;
  s1+=250.0;
 } else {
  s0+=2039.0;
  s1+=25.0;
 }
} else {
 if(i0<308.5609130859375){
  s0+=28.0;
  s1+=6.0;
 } else {
  s0+=19.0;
  s1+=1810.0;
 }
}
if(i0<399.22314453125){
 if(i1<7.672821044921875){
  s0+=2079.0;
  s1+=11.0;
 } else {
  s0+=14.0;
  s1+=71.0;
 }
} else {
 if(i0<575.8427734375){
  s0+=127.0;
  s1+=359.0;
 } else {
  s0+=33.0;
  s1+=1589.0;
 }
}
if(i1<4.690643310546875){
 if(i1<-2.376190185546875){
  s0+=2047.0;
  s1+=185.0;
 } else {
  s0+=113.0;
  s1+=71.0;
 }
} else {
 if(i1<7.429473876953125){
  s0+=19.0;
  s1+=29.0;
 } else {
  s0+=13.0;
  s1+=1806.0;
 }
}
if(i0<405.945556640625){
 if(i2<-6767.4921875){
  s0+=628.0;
  s1+=80.0;
 } else {
  s0+=1501.0;
  s1+=8.0;
 }
} else {
 if(i0<536.8392944335938){
  s0+=85.0;
  s1+=236.0;
 } else {
  s0+=29.0;
  s1+=1716.0;
 }
}
if(i2<-10154.09765625){
 if(i0<454.7654724121094){
  s0+=175.0;
  s1+=45.0;
 } else {
  s0+=102.0;
  s1+=1831.0;
 }
} else {
 if(i0<360.9933776855469){
  s0+=1954.0;
  s1+=60.0;
 } else {
  s0+=6.0;
  s1+=110.0;
 }
}
if(i0<364.408935546875){
 if(i1<7.429473876953125){
  s0+=1993.0;
  s1+=3.0;
 } else {
  s0+=8.0;
  s1+=51.0;
 }
} else {
 if(i0<609.3474731445312){
  s0+=160.0;
  s1+=452.0;
 } else {
  s0+=28.0;
  s1+=1588.0;
 }
}
if(i2<-10122.564453125){
 if(i0<435.5802917480469){
  s0+=151.0;
  s1+=16.0;
 } else {
  s0+=105.0;
  s1+=1897.0;
 }
} else {
 if(i2<-7632.1923828125){
  s0+=355.0;
  s1+=125.0;
 } else {
  s0+=1606.0;
  s1+=28.0;
 }
}
if(i1<3.3301162719726562){
 if(i2<-14161.4375){
  s0+=75.0;
  s1+=182.0;
 } else {
  s0+=2107.0;
  s1+=58.0;
 }
} else {
 if(i0<303.3450927734375){
  s0+=18.0;
  s1+=6.0;
 } else {
  s0+=16.0;
  s1+=1821.0;
 }
}
if(i1<6.4498443603515625){
 if(i0<455.9241943359375){
  s0+=2097.0;
  s1+=40.0;
 } else {
  s0+=74.0;
  s1+=267.0;
 }
} else {
 if(i0<302.26934814453125){
  s0+=10.0;
  s1+=2.0;
 } else {
  s0+=11.0;
  s1+=1782.0;
 }
}
if(i2<-11083.51953125){
 if(i0<493.781494140625){
  s0+=111.0;
  s1+=34.0;
 } else {
  s0+=70.0;
  s1+=1783.0;
 }
} else {
 if(i1<7.672821044921875){
  s0+=2036.0;
  s1+=26.0;
 } else {
  s0+=11.0;
  s1+=212.0;
 }
}
if(i0<397.90911865234375){
 if(i2<-6174.80126953125){
  s0+=725.0;
  s1+=92.0;
 } else {
  s0+=1269.0;
  s1+=1.0;
 }
} else {
 if(i2<-16973.927734375){
  s0+=19.0;
  s1+=1337.0;
 } else {
  s0+=149.0;
  s1+=691.0;
 }
}
if(i1<6.007179260253906){
 if(i2<-13803.9208984375){
  s0+=65.0;
  s1+=219.0;
 } else {
  s0+=2128.0;
  s1+=59.0;
 }
} else {
 if(i2<-7084.68310546875){
  s0+=7.0;
  s1+=1780.0;
 } else {
  s0+=16.0;
  s1+=9.0;
 }
}
if(i1<6.5342559814453125){
 if(i1<-4.3166656494140625){
  s0+=1986.0;
  s1+=187.0;
 } else {
  s0+=216.0;
  s1+=106.0;
 }
} else {
 if(i0<303.9437561035156){
  s0+=10.0;
  s1+=2.0;
 } else {
  s0+=5.0;
  s1+=1771.0;
 }
}
if(i2<-10122.12890625){
 if(i2<-13228.021484375){
  s0+=86.0;
  s1+=1573.0;
 } else {
  s0+=175.0;
  s1+=308.0;
 }
} else {
 if(i2<-7347.16162109375){
  s0+=409.0;
  s1+=133.0;
 } else {
  s0+=1582.0;
  s1+=17.0;
 }
}
if(i0<397.7173156738281){
 if(i1<7.672821044921875){
  s0+=2064.0;
  s1+=16.0;
 } else {
  s0+=17.0;
  s1+=68.0;
 }
} else {
 if(i2<-16971.7890625){
  s0+=22.0;
  s1+=1302.0;
 } else {
  s0+=130.0;
  s1+=664.0;
 }
}
if(i1<5.986244201660156){
 if(i1<-2.8140907287597656){
  s0+=2028.0;
  s1+=202.0;
 } else {
  s0+=150.0;
  s1+=95.0;
 }
} else {
 if(i0<310.71514892578125){
  s0+=17.0;
  s1+=6.0;
 } else {
  s0+=12.0;
  s1+=1773.0;
 }
}
if(i0<409.2508544921875){
 if(i2<-6877.66748046875){
  s0+=631.0;
  s1+=88.0;
 } else {
  s0+=1453.0;
  s1+=8.0;
 }
} else {
 if(i0<538.2669677734375){
  s0+=107.0;
  s1+=225.0;
 } else {
  s0+=42.0;
  s1+=1729.0;
 }
}
if(i1<3.3112564086914062){
 if(i2<-12438.15234375){
  s0+=108.0;
  s1+=213.0;
 } else {
  s0+=2050.0;
  s1+=47.0;
 }
} else {
 if(i2<-6131.056640625){
  s0+=25.0;
  s1+=1826.0;
 } else {
  s0+=14.0;
 }
}
if(i1<4.6247711181640625){
 if(i0<474.1146545410156){
  s0+=2027.0;
  s1+=25.0;
 } else {
  s0+=65.0;
  s1+=241.0;
 }
} else {
 if(i1<13.903305053710938){
  s0+=35.0;
  s1+=118.0;
 } else {
  s0+=7.0;
  s1+=1765.0;
 }
}
if(i2<-10155.67578125){
 if(i0<478.9383850097656){
  s0+=164.0;
  s1+=82.0;
 } else {
  s0+=64.0;
  s1+=1858.0;
 }
} else {
 if(i0<356.0542907714844){
  s0+=1921.0;
  s1+=45.0;
 } else {
  s0+=20.0;
  s1+=129.0;
 }
}
if(i1<3.3112564086914062){
 if(i0<474.216064453125){
  s0+=2094.0;
  s1+=43.0;
 } else {
  s0+=69.0;
  s1+=202.0;
 }
} else {
 if(i2<-6131.056640625){
  s0+=25.0;
  s1+=1833.0;
 } else {
  s0+=17.0;
 }
}
if(i2<-10154.09765625){
 if(i1<-18.234619140625){
  s0+=214.0;
  s1+=89.0;
 } else {
  s0+=23.0;
  s1+=1773.0;
 }
} else {
 if(i2<-7278.2216796875){
  s0+=439.0;
  s1+=175.0;
 } else {
  s0+=1553.0;
  s1+=17.0;
 }
}
if(i2<-11418.033203125){
 if(i1<-44.92584228515625){
  s0+=98.0;
  s1+=33.0;
 } else {
  s0+=42.0;
  s1+=1708.0;
 }
} else {
 if(i0<360.9385986328125){
  s0+=2043.0;
  s1+=46.0;
 } else {
  s0+=71.0;
  s1+=242.0;
 }
}
if(i0<397.752197265625){
 if(i1<8.522415161132812){
  s0+=2103.0;
  s1+=15.0;
 } else {
  s0+=14.0;
  s1+=61.0;
 }
} else {
 if(i2<-18610.484375){
  s0+=14.0;
  s1+=1124.0;
 } else {
  s0+=128.0;
  s1+=824.0;
 }
}
if(i1<3.3112564086914062){
 if(i1<-4.3166656494140625){
  s0+=1958.0;
  s1+=167.0;
 } else {
  s0+=175.0;
  s1+=74.0;
 }
} else {
 if(i0<310.8768310546875){
  s0+=32.0;
  s1+=11.0;
 } else {
  s0+=17.0;
  s1+=1849.0;
 }
}
if(i0<401.87628173828125){
 if(i1<6.9555206298828125){
  s0+=2108.0;
  s1+=9.0;
 } else {
  s0+=11.0;
  s1+=73.0;
 }
} else {
 if(i0<454.5257568359375){
  s0+=51.0;
  s1+=88.0;
 } else {
  s0+=84.0;
  s1+=1859.0;
 }
}
if(i0<404.2063903808594){
 if(i0<320.4058837890625){
  s0+=1866.0;
  s1+=15.0;
 } else {
  s0+=176.0;
  s1+=92.0;
 }
} else {
 if(i1<-23.91558837890625){
  s0+=110.0;
  s1+=85.0;
 } else {
  s0+=27.0;
  s1+=1912.0;
 }
}
if(i1<3.3112564086914062){
 if(i2<-13803.9208984375){
  s0+=66.0;
  s1+=184.0;
 } else {
  s0+=2069.0;
  s1+=45.0;
 }
} else {
 if(i1<7.429473876953125){
  s0+=26.0;
  s1+=62.0;
 } else {
  s0+=18.0;
  s1+=1813.0;
 }
}
if(i0<397.90911865234375){
 if(i1<6.007179260253906){
  s0+=2004.0;
  s1+=15.0;
 } else {
  s0+=17.0;
  s1+=85.0;
 }
} else {
 if(i1<-18.281112670898438){
  s0+=107.0;
  s1+=110.0;
 } else {
  s0+=15.0;
  s1+=1930.0;
 }
}
if(i2<-10154.09765625){
 if(i0<454.7654724121094){
  s0+=160.0;
  s1+=34.0;
 } else {
  s0+=95.0;
  s1+=1872.0;
 }
} else {
 if(i2<-7611.494140625){
  s0+=320.0;
  s1+=148.0;
 } else {
  s0+=1630.0;
  s1+=24.0;
 }
}
if(i0<367.86083984375){
 if(i0<321.56256103515625){
  s0+=1892.0;
  s1+=6.0;
 } else {
  s0+=131.0;
  s1+=48.0;
 }
} else {
 if(i2<-13367.9609375){
  s0+=58.0;
  s1+=1646.0;
 } else {
  s0+=108.0;
  s1+=394.0;
 }
}
if(i0<406.1541442871094){
 if(i2<-6763.92822265625){
  s0+=633.0;
  s1+=87.0;
 } else {
  s0+=1439.0;
  s1+=10.0;
 }
} else {
 if(i0<608.3754272460938){
  s0+=119.0;
  s1+=375.0;
 } else {
  s0+=28.0;
  s1+=1592.0;
 }
}
if(i1<6.007179260253906){
 if(i2<-12552.05859375){
  s0+=112.0;
  s1+=205.0;
 } else {
  s0+=2018.0;
  s1+=53.0;
 }
} else {
 if(i0<310.71514892578125){
  s0+=11.0;
  s1+=7.0;
 } else {
  s0+=15.0;
  s1+=1862.0;
 }
}
if(i2<-10416.853515625){
 if(i2<-12410.908203125){
  s0+=104.0;
  s1+=1639.0;
 } else {
  s0+=122.0;
  s1+=170.0;
 }
} else {
 if(i1<7.429473876953125){
  s0+=2018.0;
  s1+=18.0;
 } else {
  s0+=10.0;
  s1+=202.0;
 }
}
if(i2<-9481.560546875){
 if(i1<-18.176315307617188){
  s0+=282.0;
  s1+=118.0;
 } else {
  s0+=16.0;
  s1+=1814.0;
 }
} else {
 if(i2<-6877.66748046875){
  s0+=459.0;
  s1+=130.0;
 } else {
  s0+=1456.0;
  s1+=8.0;
 }
}
if(i1<2.8080406188964844){
 if(i1<-1.9239578247070312){
  s0+=2041.0;
  s1+=200.0;
 } else {
  s0+=108.0;
  s1+=51.0;
 }
} else {
 if(i2<-7599.6337890625){
  s0+=13.0;
  s1+=1803.0;
 } else {
  s0+=36.0;
  s1+=31.0;
 }
}
if(i2<-10154.09765625){
 if(i2<-12908.716796875){
  s0+=93.0;
  s1+=1642.0;
 } else {
  s0+=168.0;
  s1+=234.0;
 }
} else {
 if(i2<-8244.673828125){
  s0+=250.0;
  s1+=104.0;
 } else {
  s0+=1742.0;
  s1+=50.0;
 }
}
if(i1<5.3414154052734375){
 if(i1<-1.6529617309570312){
  s0+=2038.0;
  s1+=216.0;
 } else {
  s0+=86.0;
  s1+=69.0;
 }
} else {
 if(i1<15.111175537109375){
  s0+=25.0;
  s1+=126.0;
 } else {
  s0+=2.0;
  s1+=1721.0;
 }
}
if(i1<5.154563903808594){
 if(i0<492.66339111328125){
  s0+=2123.0;
  s1+=42.0;
 } else {
  s0+=54.0;
  s1+=213.0;
 }
} else {
 if(i0<308.0574645996094){
  s0+=15.0;
  s1+=5.0;
 } else {
  s0+=12.0;
  s1+=1819.0;
 }
}
if(i1<6.4498443603515625){
 if(i0<493.47802734375){
  s0+=2149.0;
  s1+=30.0;
 } else {
  s0+=61.0;
  s1+=239.0;
 }
} else {
 if(i2<-6503.6640625){
  s0+=12.0;
  s1+=1783.0;
 } else {
  s0+=4.0;
  s1+=5.0;
 }
}
if(i1<5.02783203125){
 if(i1<-4.8751220703125){
  s0+=1926.0;
  s1+=165.0;
 } else {
  s0+=180.0;
  s1+=86.0;
 }
} else {
 if(i0<280.22613525390625){
  s0+=8.0;
  s1+=2.0;
 } else {
  s0+=14.0;
  s1+=1902.0;
 }
}
 float max_s=s0;
 int cls=1;
 if (max_s < s1) {
  max_s = s1;
  cls=2;
 }
 WRITE_IMAGE (out, POS_out_INSTANCE(x,y,z,0), cls);
}
