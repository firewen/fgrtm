gmt begin seis png
	gmt basemap -JX10c/5c -R0/122.85/0/4 -BWSen -Bxa20+l"Time(s)" -By+l"St. No."
	gmt sac -M1c -S30c -T+t-5 << EOF
GZ.BJT.R.sac	0 1.0
GZ.BJT.T.sac	0 2.0
GZ.BJT.Z.sac	0 3.0
EOF
gmt end

