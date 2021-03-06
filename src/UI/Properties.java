package UI;/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Kiky
 */
import semiglobalalignment.SemiGlobalAlignment;
public class Properties extends javax.swing.JFrame {
    /**
     * Creates new form Properties
     */
    public Properties() {
        initComponents();
	errorGap.setVisible(false);
	errorGapType.setVisible(false);
	errorGapInisialisasi.setVisible(false);
	errorGapIterasi.setVisible(false);
	errorMatriks.setVisible(false);
	//public static int pilihangap, pilihanmatriks;
    }

  /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        bindingGroup = new org.jdesktop.beansbinding.BindingGroup();

        GapLinear = new javax.swing.JRadioButton();
        GapAffine = new javax.swing.JRadioButton();
        LabelGap = new javax.swing.JLabel();
        PAM250 = new javax.swing.JRadioButton();
        BLOSUM62 = new javax.swing.JRadioButton();
        Labelmatriks = new javax.swing.JLabel();
        OK = new javax.swing.JButton();
        Gap = new javax.swing.JTextField();
        GapInisialisasi = new javax.swing.JTextField();
        GapIterasi = new javax.swing.JTextField();
        cancel = new javax.swing.JButton();
        errorGap = new javax.swing.JLabel();
        errorGapInisialisasi = new javax.swing.JLabel();
        errorGapIterasi = new javax.swing.JLabel();
        errorMatriks = new javax.swing.JLabel();
        errorGapType = new javax.swing.JLabel();
        bQ = new javax.swing.JCheckBox();
        eQ = new javax.swing.JCheckBox();
        bD = new javax.swing.JCheckBox();
        eD = new javax.swing.JCheckBox();
        jLabel1 = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Properties");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        GapLinear.setSelected(false);
        GapLinear.setText("Linear Gap ");

        org.jdesktop.beansbinding.Binding binding = org.jdesktop.beansbinding.Bindings.createAutoBinding(org.jdesktop.beansbinding.AutoBinding.UpdateStrategy.READ_WRITE, Gap, org.jdesktop.beansbinding.ELProperty.create("${enabled}"), GapLinear, org.jdesktop.beansbinding.BeanProperty.create("actionCommand"));
        bindingGroup.addBinding(binding);

        GapLinear.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                GapLinearMouseClicked(evt);
            }
        });

        GapAffine.setText("Affine Gap");
        GapAffine.setToolTipText("");
        GapAffine.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                GapAffineMouseClicked(evt);
            }
        });
        GapAffine.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                GapLinear(evt);
            }
        });

        LabelGap.setFont(new java.awt.Font("Tahoma", 1, 12)); // NOI18N
        LabelGap.setText("Gap  ");

        PAM250.setText("PAM250");
        PAM250.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                PAM250MouseClicked(evt);
            }
        });

        BLOSUM62.setText("BLOSUM62");
        BLOSUM62.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                BLOSUM62MouseClicked(evt);
            }
        });

        Labelmatriks.setFont(new java.awt.Font("Tahoma", 1, 12)); // NOI18N
        Labelmatriks.setText("Matrix ");

        OK.setFont(new java.awt.Font("Tahoma", 1, 11)); // NOI18N
        OK.setText("OK");
        OK.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                OKMouseClicked(evt);
            }
        });

        Gap.setText("Gap");
        Gap.setEnabled(GapLinear.isSelected());
        Gap.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                GapMouseClicked(evt);
            }
        });
        Gap.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                GapFocusLost(evt);
            }
        });

        GapInisialisasi.setText("Initialize gap");
        GapInisialisasi.setEnabled(GapAffine.isSelected());
        GapInisialisasi.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                GapInisialisasiMouseClicked(evt);
            }
        });
        GapInisialisasi.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                GapInisialisasiFocusLost(evt);
            }
        });

        GapIterasi.setText("Iteration gap");
        GapIterasi.setEnabled(GapAffine.isSelected());
        GapIterasi.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                GapIterasiMouseClicked(evt);
            }
        });
        GapIterasi.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                GapIterasiFocusLost(evt);
            }
        });

        cancel.setText("Cancel");
        cancel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                cancelMouseClicked(evt);
            }
        });

        errorGap.setForeground(new java.awt.Color(204, 0, 0));
        errorGap.setLabelFor(Gap);
        errorGap.setText("this field cannot be empty");

        errorGapInisialisasi.setForeground(new java.awt.Color(204, 0, 0));
        errorGapInisialisasi.setText("this field cannot be empty");

        errorGapIterasi.setForeground(new java.awt.Color(204, 0, 0));
        errorGapIterasi.setText("this field cannot be empty");

        errorMatriks.setForeground(new java.awt.Color(204, 0, 0));
        errorMatriks.setText("please choose a matrix");

        errorGapType.setForeground(new java.awt.Color(204, 0, 0));
        errorGapType.setText("please choose a gap type");

        bQ.setText("beginning of sequence q");
        bQ.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                bQMouseClicked(evt);
            }
        });

        eQ.setText("end of sequence q");

        bD.setText("beginning of sequence d");
        bD.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bDActionPerformed(evt);
            }
        });

        eD.setText("end of sequence d");
        eD.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eDActionPerformed(evt);
            }
        });

        jLabel1.setFont(new java.awt.Font("Tahoma", 1, 11)); // NOI18N
        jLabel1.setText("Ungraded Gap");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(27, 27, 27)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(22, 22, 22)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(errorGapInisialisasi, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addGap(225, 225, 225))
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(layout.createSequentialGroup()
                                                .addComponent(GapInisialisasi, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addGap(81, 81, 81)
                                                .addComponent(GapIterasi, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE))
                                            .addGroup(layout.createSequentialGroup()
                                                .addGap(106, 106, 106)
                                                .addComponent(BLOSUM62))
                                            .addGroup(layout.createSequentialGroup()
                                                .addGap(171, 171, 171)
                                                .addComponent(errorGapIterasi))
                                            .addGroup(layout.createSequentialGroup()
                                                .addComponent(Gap, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                                .addComponent(errorGap)))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(LabelGap)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(errorGapType))
                                    .addComponent(GapAffine)
                                    .addComponent(GapLinear))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(OK, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)))
                .addComponent(cancel, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(28, 28, 28))
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(29, 29, 29)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(PAM250)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(Labelmatriks)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(errorMatriks))))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(418, 418, 418)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel1)
                            .addComponent(eQ)
                            .addComponent(eD)
                            .addComponent(bQ)
                            .addComponent(bD))))
                .addContainerGap(49, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(26, 26, 26)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(LabelGap)
                    .addComponent(errorGapType))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(GapLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(Gap, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(errorGap))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GapAffine)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(GapIterasi, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(GapInisialisasi, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(errorGapInisialisasi, javax.swing.GroupLayout.PREFERRED_SIZE, 14, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(errorGapIterasi)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(bQ)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eQ)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(bD)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eD)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED, 15, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(errorMatriks)
                    .addComponent(Labelmatriks))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(BLOSUM62)
                    .addComponent(PAM250))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cancel)
                    .addComponent(OK))
                .addGap(35, 35, 35))
        );

        Gap.getAccessibleContext().setAccessibleName("");

        bindingGroup.bind();

        setBounds(100, 100, 626, 315);
    }// </editor-fold>//GEN-END:initComponents

    private void GapLinear(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_GapLinear
        // TODO add your handling code here:
    }//GEN-LAST:event_GapLinear

    private void GapLinearMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_GapLinearMouseClicked
        errorGapType.setVisible(false);errorGap.setVisible(false);
	errorGapInisialisasi.setVisible(false);errorGapIterasi.setVisible(false);
	if (GapLinear.isSelected()==true){
	    Gap.setEnabled(true);
	    GapAffine.setSelected(false);
	    GapIterasi.setEnabled(false);
	    GapIterasi.setText("Iteration gap");
	    GapInisialisasi.setEnabled(false);
	    GapInisialisasi.setText("Initialize gap");
	}
	else if (GapLinear.isSelected()==false){
	    Gap.setEnabled(false);
	    GapAffine.setSelected(true);
	    GapIterasi.setEnabled(true);
	    GapInisialisasi.setEnabled(true);
	    Gap.setText("Gap");
	}
    }//GEN-LAST:event_GapLinearMouseClicked

    private void cancelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_cancelMouseClicked
        semiglobalalignment.SemiGlobalAlignment.fungsi.prim.setEnabled(true);
	this.dispose();// TODO add your handling code here:
    }//GEN-LAST:event_cancelMouseClicked

    private void GapAffineMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_GapAffineMouseClicked
        errorGapType.setVisible(false);errorGapInisialisasi.setVisible(false);
	errorGapIterasi.setVisible(false);errorGap.setVisible(false);
	if (GapAffine.isSelected()==true){
	    GapIterasi.setEnabled(true);
	    GapInisialisasi.setEnabled(true);
	    GapLinear.setSelected(false);
	    Gap.setEnabled(false);
	    Gap.setText("Gap");
	}
	else if (GapAffine.isSelected()==false){
	    GapIterasi.setEnabled(false);
	    GapInisialisasi.setEnabled(false);
	    GapLinear.setSelected(true);
	    Gap.setEnabled(true);
	    GapIterasi.setText("Iteration gap");
	    GapInisialisasi.setText("Initialize gap");
	}
	// TODO add your handling code here:
    }//GEN-LAST:event_GapAffineMouseClicked

    private void PAM250MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_PAM250MouseClicked
        errorMatriks.setVisible(false);
	if (PAM250.isSelected()==true)
	    BLOSUM62.setSelected(false);
	else
	    BLOSUM62.setSelected(true);// TODO add your handling code here:
    }//GEN-LAST:event_PAM250MouseClicked

    private void BLOSUM62MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_BLOSUM62MouseClicked
        errorMatriks.setVisible(false); 
	if (BLOSUM62.isSelected()==true)
	    PAM250.setSelected(false);
	 else
	     PAM250.setSelected(true);// TODO add your handling code here:
    }//GEN-LAST:event_BLOSUM62MouseClicked

    private void GapMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_GapMouseClicked
        errorGap.setVisible(false);
	if (Gap.isFocusOwner()==true && GapLinear.isSelected()==true && Gap.getText().equals("Gap"))
	    Gap.setText(null);
	// TODO add your handling code here:
    }//GEN-LAST:event_GapMouseClicked

    private void GapInisialisasiMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_GapInisialisasiMouseClicked
        errorGapInisialisasi.setVisible(false);
	if (GapInisialisasi.isFocusOwner()==true && GapAffine.isSelected()==true && GapInisialisasi.getText().equals("Initialize gap"))
	    GapInisialisasi.setText(null);// TODO add your handling code here:
    }//GEN-LAST:event_GapInisialisasiMouseClicked

    private void GapIterasiMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_GapIterasiMouseClicked
        errorGapIterasi.setVisible(false);
	if (GapIterasi.isFocusOwner()==true && GapAffine.isSelected()==true && GapIterasi.getText().equals("Iteration gap"))
	    GapIterasi.setText(null);
	
    }//GEN-LAST:event_GapIterasiMouseClicked

    private void GapIterasiFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GapIterasiFocusLost
        errorGapIterasi.setVisible(false);
	if (GapIterasi.getText().equals(""))
	    GapIterasi.setText("Iteration gap");// TODO add your handling code here:
    }//GEN-LAST:event_GapIterasiFocusLost

    private void GapInisialisasiFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GapInisialisasiFocusLost
        if (GapInisialisasi.getText().equals(""))
	    GapInisialisasi.setText("Initialize gap");// TODO add your handling code here:
    }//GEN-LAST:event_GapInisialisasiFocusLost

    private void GapFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GapFocusLost
        if (Gap.getText().equals(""))
	    Gap.setText("Gap");// TODO add your handling code here:
    }//GEN-LAST:event_GapFocusLost

    private void OKMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_OKMouseClicked
	if (GapLinear.isSelected()==false && GapAffine.isSelected()==false)
	    errorGapType.setVisible(true);
	if (GapLinear.isSelected()==true && (Gap.getText().equals("Gap")||Gap.getText().equals("")))
	    errorGap.setVisible(true);
	if (GapAffine.isSelected()==true){
	    if (GapIterasi.getText().equals("")||GapIterasi.getText().equals("Iteration gap"))
		errorGapIterasi.setVisible(true);
	    if (GapInisialisasi.getText().equals("")||GapInisialisasi.getText().equals("Initialize gap"))
		errorGapInisialisasi.setVisible(true);
	}
	if (PAM250.isSelected()==false && BLOSUM62.isSelected()==false)
	    errorMatriks.setVisible(true);
	if (errorGap.isVisible()==false&&errorGapType.isVisible()==false&&errorGapInisialisasi.isVisible()==false&&errorMatriks.isVisible()==false){
	    semiglobalalignment.SemiGlobalAlignment.fungsi.inisialisasi();
	    if (GapLinear.isSelected()==true){
		semiglobalalignment.SemiGlobalAlignment.fungsi.pilihangap=1;
		semiglobalalignment.SemiGlobalAlignment.fungsi.gapIterasi=Integer.parseInt(Gap.getText());
		semiglobalalignment.SemiGlobalAlignment.fungsi.gapInisialisasi=semiglobalalignment.SemiGlobalAlignment.fungsi.gapIterasi;
	    }
	    else{
		semiglobalalignment.SemiGlobalAlignment.fungsi.pilihangap=2;
		semiglobalalignment.SemiGlobalAlignment.fungsi.gapIterasi=Integer.parseInt(GapIterasi.getText());
		semiglobalalignment.SemiGlobalAlignment.fungsi.gapInisialisasi=Integer.parseInt(GapInisialisasi.getText());
	    }
	    if (PAM250.isSelected()==true)
		semiglobalalignment.SemiGlobalAlignment.fungsi.pilihanmatriks=1;
	    else
		semiglobalalignment.SemiGlobalAlignment.fungsi.pilihanmatriks=2;
	    if (bQ.isSelected()==true)
		SemiGlobalAlignment.fungsi.awalQ=true;
	    if (bD.isSelected()==true)
		SemiGlobalAlignment.fungsi.awalD=true;
	    if (eQ.isSelected()==true)
		SemiGlobalAlignment.fungsi.akhirQ=true;
	    if (eD.isSelected()==true)
		SemiGlobalAlignment.fungsi.akhirD=true;
	    semiglobalalignment.SemiGlobalAlignment.fungsi.prim.setEnabled(true);
	    this.dispose();
	}
	// TODO add your handling code here:
    }//GEN-LAST:event_OKMouseClicked

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        semiglobalalignment.SemiGlobalAlignment.fungsi.prim.setEnabled(true);// TODO add your handling code here:
    }//GEN-LAST:event_formWindowClosing

    private void bDActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bDActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_bDActionPerformed

    private void eDActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eDActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_eDActionPerformed

    private void bQMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_bQMouseClicked
        if (bQ.isSelected()==true)
	    SemiGlobalAlignment.fungsi.awalQ=true;// TODO add your handling code here:
    }//GEN-LAST:event_bQMouseClicked

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(Properties.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(Properties.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(Properties.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Properties.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
	    @Override
            public void run() {
                //new Properties().setVisible(true);		
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JRadioButton BLOSUM62;
    private javax.swing.JTextField Gap;
    private javax.swing.JRadioButton GapAffine;
    private javax.swing.JTextField GapInisialisasi;
    private javax.swing.JTextField GapIterasi;
    private javax.swing.JRadioButton GapLinear;
    private javax.swing.JLabel LabelGap;
    private javax.swing.JLabel Labelmatriks;
    private javax.swing.JButton OK;
    private javax.swing.JRadioButton PAM250;
    private javax.swing.JCheckBox bD;
    private javax.swing.JCheckBox bQ;
    private javax.swing.JButton cancel;
    private javax.swing.JCheckBox eD;
    private javax.swing.JCheckBox eQ;
    private static javax.swing.JLabel errorGap;
    private static javax.swing.JLabel errorGapInisialisasi;
    private static javax.swing.JLabel errorGapIterasi;
    private static javax.swing.JLabel errorGapType;
    private static javax.swing.JLabel errorMatriks;
    private javax.swing.JLabel jLabel1;
    private org.jdesktop.beansbinding.BindingGroup bindingGroup;
    // End of variables declaration//GEN-END:variables
}
