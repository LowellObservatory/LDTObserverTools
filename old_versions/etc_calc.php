<?php
$vMode = $_POST["mode"];
$vA_Mag = $_POST["a_mag"];
$vA_SN = $_POST["a_sn"];
$vB_Mag = $_POST["b_mag"];
$vB_Exp = $_POST["b_exp"];
$vC_SN = $_POST["c_sn"];
$vC_Exp = $_POST["c_exp"];
$vBand = $_POST["band"];
$vTel = $_POST["tel"];
$vBin = $_POST["bin"];
$vSeeing = $_POST["seeing"];
$vPhase = $_POST["phase"];
$vAirmass = $_POST["airmass"];
if (!isset($_POST['submit'])) { 
?>

<html>
<head>
  <title>Lowell Observatory Exposure Time Calculator (ETC)</title>
</head>

<body style="width: 800px;">
  <h1>Lowell Observatory Exposure Time Calculator (ETC)</h1>

  <p><b>Input any two (S/N, magnitude, exposure time) and obtain the remaining value.</b>


  <form method="post" action="<?php echo $PHP_SELF;?>">
  
  	<p><b>Mode:</b><br>
  	
  	<input type="radio" name="mode" value="exptime" checked> Exposure time for a given S/N and magnitude. Please provide:<br> 
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	Magnitude: <input type="text" name="a_mag"><br>
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	S/N ratio (> 0.1): <input type="text" name="a_sn"><br><br>
  	
    <input type="radio" name="mode" value="sn"> S/N for a given exposure time and magnitude. Please provide:<br>
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	Magnitude: <input type="text" name="b_mag"><br>
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	Exposure Time (s) (> 0.001): <input type="text" name="b_exp"><br><br>
  	
    <input type="radio" name="mode" value="mag"> Magnitude for a given S/N and exposure time. Please provide:<br>
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	S/N ratio (> 0.1): <input type="text" name="c_sn"><br>
  	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	Exposure Time (s) (> 0.001): <input type="text" name="c_exp"><br>
  	
  	<p><b>Bandpass:</b>
  	<select name="band">
  		<option>U</option>
  		<option>B</option>
  		<option selected>V</option>  		
  		<option>R</option>
  		<option>I</option>
        <option>VR</option>
        <option>H&#945;-on</option>
        <option>H&#945;-off</option>
        <option>[OIII]</option>
        <option>WR</option>
        <option>u&#039;</option>
        <option>g&#039;</option>
        <option>r&#039;</option>
        <option>i&#039;</option>
        <option>z&#039;</option>
        <option>yish</option>
  	</select>  
    
	<p><b>Telescope/Instrument:</b><br>
	<input type="radio" name="tel" value="DCT" checked> DCT/LMI<br>
	<input type="radio" name="tel" value="Hall"> Hall/NASA42<br>
  
  	<p><b>CCD Binning Factor:</b>
  	<select name="bin">
  		<option>1</option>
  		<option selected>2</option>
  		<option>3</option>
  		<option>4</option>  		  		  
  	</select>
  	
  	<p><b>Seeing (arcseconds):</b>
  	<select name="seeing">
  		<option>0.6</option>
  		<option>0.8</option>
  		<option selected>1.0</option>
  		<option>1.2</option>
  		<option>1.5</option>
  		<option>2.0</option>
  		<option>2.5</option>
  	</select>
  	
  	<p><b>Lunar Phase (days):</b>
	<select name="phase">
		<option value="0">0 (new moon)</option>
		<option value="1">1</option>
		<option value="2">2</option>
		<option value="3">3</option>
		<option value="4">4</option>
		<option value="5">5</option>
		<option value="6">6</option>
		<option value="7">7 (first/last quarter)</option>
		<option value="8">8</option>
		<option value="9">9</option>
		<option value="10">10</option>
		<option value="11">11</option>
		<option value="12">12</option>
		<option value="13">13</option>
		<option value="14">14 (full moon)</option>
	</select>						
				
	<p><b>Airmass:</b>
	<select name="airmass">
		<option>1.0</option>
		<option>1.2</option>
		<option>1.4</option>
		<option>1.6</option>
		<option>1.8</option>
		<option>2.0</option>
		<option>2.2</option>
		<option>2.4</option>														
  	</select>		  		

    <p>
<input type="submit" value="submit" name="submit"><br>
<input type="reset"><br><br>
<p> Created by Phil Massey and Kathryn Neugent</p>
  <p><b><a href='ETCmethod.pdf' target="_blank">Basic equations and assumptions can be found here.</a>
   </form>
   <?
   } else {
   if ($vTel=="DCT") {
   		$scale=0.12;
   		$read_noise=6.0;
//Count ratess in e/sec/image at X=0 for a 20th mag star measured with a radius=1.4xfwhm in pixels
   		if ($vBand=="U") $Star20=120;
   		if ($vBand=="B") $Star20=770;
   		if ($vBand=="V") $Star20=670;
   		if ($vBand=="R") $Star20=630;
                if ($vBand=="I") $Star20=430;
                if ($vBand=="VR") $Star20=990;
                if ($vBand=="H&#945;-on") $Star20=15;
                if ($vBand=="H&#945;-off")$Star20=70;
                if ($vBand=="[OIII]") $Star20=18;
                if ($vBand=="WR") $Star20=35;
                if ($vBand=="u'") $Star20=220;
                if ($vBand=="g'") $Star20=900;
                if ($vBand=="r'") $Star20=750;
                if ($vBand=="i'") $Star20=520;
                if ($vBand=="z'") $Star20=325;
                if ($vBand=="yish") $Star20=20;
   } else {
   		$scale=0.37;
   		$read_noise=4.6;
   		if ($vBand=="U") $Star20=3.1;
   		if ($vBand=="B") $Star20=17.1;
   		if ($vBand=="V") $Star20=17.5;
   		if ($vBand=="R") $Star20=19.7;
   		if ($vBand=="I") $Star20=12.8;
                if ($vBand=="VR") $Star20=25.0;
                if ($vBand=="H&#945;-on") $Star20=0.4;
                if ($vBand=="H&#945;-off") $Star20=2.0;
   }

   // Determine sky counts per pixel per second:
   
   if ($vBand=="U"|$vBand=="u'") {
   		$extinction=0.55;
   		$sky0=22.0;
   		$sky1=-0.2666;
   		$sky2=-0.0760;
   } elseif ($vBand=="B" | $vBand=="g'") {
   		$extinction=0.25;
   		$sky0=22.7;
   		$sky1=-0.0998;
   		$sky2=-0.00953;
		$extinction=0.20;
   } elseif ($vBand=="WR") {
                $extinction=0.20;
                $sky0=22.5;
                $sky1=-0.0998;
                $sky2=-0.00953;
   } elseif ($vBand=="[OIII]") {
                $extinction=0.20;
                $sky0=22.3;
                $sky1=-0.00998;
                $sky2=-0.00953;
   } elseif ($vBand=="V") {
   		$extinction=0.14;
   		$sky0=21.8;
   		$sky1=-0.0153;
   		$sky2=-0.00838;
   } elseif ($vBand=="R" | $vBand=="r'") {
   		$extinction=0.10;
   		$sky0=20.9;
   		$sky1=-0.0211;
   		$sky2=-0.00364;
   } elseif ($vBand=="I"|$vBand=="i'"|$vBand=="yish"|$vBand=="z'") {
   		$extinction=0.05;
   		$sky0=19.9;
   		$sky1=-0.018;
   		$sky2=-0.005;
   } elseif ($vBand=="VR") {
       $extinction=0.12;
       $sky0=21.4;
       $sky1=-0.011;
       $sky2=-0.00311;
   } elseif ($vBand=="H&#945;-on") {
   		$extinction=0.10;
   		$sky0=20.9;
   		$sky1=-0.0211;
   		$sky2=-0.00364;
   } elseif ($vBand=="H&#945;-off") {
   		$extinction=0.10;
   		$sky0=20.9;
   		$sky1=-0.0211;
   		$sky2=-0.00364;
   }
   $sky_brightness_per_arcsec2=$sky0+$sky1*$vPhase+$sky2*pow($vPhase,2);
   $sky_count_per_arcsec2_per_sec=$Star20*pow(10,-(($sky_brightness_per_arcsec2-20)/2.5));
       $rscale=$scale*$vBin;
       $sky_count_per_pixel_per_sec=$sky_count_per_arcsec2_per_sec*pow($rscale,2);
       $fwhm=$vSeeing/$rscale;
   $Number_pixels=1.4*pow($fwhm,2);
       if($Number_pixels < 9) $Number_pixels=9;
   $sky_count_per_sec_per_ap=$Number_pixels * $sky_count_per_pixel_per_sec;

   // Calculate read-noise contribution:
   
       $read_contribution=$read_noise*sqrt($Number_pixels);

   // Now calculate stuff depending upon the mode:
   if ($vMode=="exptime") {
   		$mag_corrected=$vA_Mag+$extinction*($vAirmass);
   		$counts_from_star_per_sec=$Star20*pow(10,-(($mag_corrected-20)/2.5));
        $sn=$vA_SN;
        $sn2=pow($vA_SN,2);
   		$A=pow($counts_from_star_per_sec,2);
   		$B=-$sn2*($counts_from_star_per_sec+$sky_count_per_sec_per_ap);
        $C=-$sn2*pow($read_contribution,2);
   		$exptime=(-$B+sqrt($B*$B-4.0*$A*$C))/(2.0*$A);
   		$mag=$vA_Mag;
   } elseif ($vMode=="sn") {
   		$mag_corrected=$vB_Mag+$extinction*($vAirmass);
   		$counts_from_star_per_sec=$Star20*pow(10,-(($mag_corrected-20)/2.5));
   		$SNR=$counts_from_star_per_sec*$vB_Exp/ (sqrt($counts_from_star_per_sec*$vB_Exp+$sky_count_per_sec_per_ap*$vB_Exp+pow($read_contribution,2)));
   		$mag=$vB_Mag;
   		$sn=$SNR;
   		$exptime=$vB_Exp;
   } elseif ($vMode=="mag") {
   		$A=pow($vC_Exp,2);
   		$B=-$vC_Exp*pow($vC_SN,2);
   		$C=-pow($vC_SN,2)*($sky_count_per_sec_per_ap*$vC_Exp+pow($read_contribution,2));
   		$counts_from_star_per_sec=(-$B+sqrt(pow($B,2)-4*$A*$C))/(2*$A);
   		$mag_raw=-2.5*log10($counts_from_star_per_sec/$Star20)+20.0;
   		$mag=$mag_raw-$extinction*($vAirmass);
   		$sn=$vC_SN;
   		$exptime=$vC_Exp;
   }  
       $peak=$counts_from_star_per_sec*$exptime/(1.13*pow($fwhm,2));
echo "<html>";
echo "<head>";
echo "<title>Exposure Time Calculator (ETC)</title>";
echo "</head>";
echo "<h1>Exposure Time Calculator (ETC) Results</h1>";

echo "<b>Bandpass =</b> $vBand<br>";
$mag = round($mag,2);
echo "<b>Magnitude =</b> $mag<br>";
       $nc_sn=round($sn,1);
echo "<b>S/N =</b> $nc_sn<br>";
$nc_exptime=round($exptime,3);
echo "<b>Exposure Time (s) =</b> $nc_exptime<br>";
echo "<b>Binning =</b> $vBin<br><br>";

$Number_pixel = round($Number_pixels,2);
echo "<b>Number of pixels in measuring aperture =</b> $Number_pixel<br><br>";
$sky_count_per_pixel = round($sky_count_per_pixel_per_sec*$exptime,1);
echo "<b>Sky brightness (e per pixel) =</b> $sky_count_per_pixel <br><br><br>";
       $Num_star=round($counts_from_star_per_sec*$exptime,1);
echo "<b>Number of e from star=</b> $Num_star<br>";
       $peak=round($peak,0);
echo "<b>Peak number of e from star=</b> $peak<br><br>";

$nc_star=round(sqrt($counts_from_star_per_sec*$exptime),2);
echo "<b>Noise contribution (Star) [e] =</b> $nc_star<br>";
$nc_sky=round(sqrt($sky_count_per_sec_per_ap*$exptime),2);
echo "<b>Noise contribution (Sky) [e] =</b> $nc_sky<br>";
$read_contribution=round($read_contribution,2);
echo "<b>Noise contribution (CCD) [e] =</b> $read_contribution<br>";
}
?>  
